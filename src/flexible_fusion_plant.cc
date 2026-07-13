#include "flexible_fusion_plant.h"

using cyclus::CompMap;
using cyclus::Composition;
using cyclus::DoubleDistribution;
using cyclus::FixedDoubleDist;
using cyclus::FixedIntDist;
using cyclus::IntDistribution;
using cyclus::KeyError;
using cyclus::Material;
using cyclus::toolkit::ResBuf;

namespace tricycle {

const int tritium_id = 10030000;
const cyclus::CompMap T = {{tritium_id, 1}};
const cyclus::Composition::Ptr tritium_comp = cyclus::Composition::CreateFromAtom(T);
const double mass_tritium = pyne::atomic_mass(tritium_id) / (1000.0 * pyne::N_A);

const double MW_to_W = 1000000;
const double MeV_to_J = 1.6021766E-13;
const double energy_DT = 17.6 * MeV_to_J;

// Defined to avoid using eps_rsrc(): this has a value of 1e-6.
// That would risk missing/losing a significant amount of tritium
// which can be tracked in micrograms.
const double tritium_eps = 1.0e-12;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
FlexibleFusionPlant::FlexibleFusionPlant(cyclus::Context* ctx)
    : cyclus::Facility(ctx) {

  tritium_storage = ResBuf<Material>(true);
  tritium_excess = ResBuf<Material>(true);

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::string FlexibleFusionPlant::str() {
  return Facility::str();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void FlexibleFusionPlant::EnterNotify() {
  cyclus::Facility::EnterNotify();

  // Only used to instantiate the tracker
  double fuel_limit = 1000.0;
  fuel_tracker.Init({&tritium_storage}, fuel_limit);

  sequestered_tritium = cyclus::Material::CreateUntracked(0.0, tritium_comp);
  
  burn_rate = mass_tritium * fusion_power * MW_to_W/ 
	  (conversion_efficiency * energy_DT);
  feed_rate = burn_rate / TBE;

  failure_probability = 1.0 - std::exp(-failure_frequency * context()->dt() / cyclusYear);  
  
  int N = components.size();
  
  // Size vector of intermediate buffers
  tritium_elsewhere = std::vector<ResBuf<Material>>(N, ResBuf<Material>(true));

  // Build component lookup map
  for (int i = 0; i < N; ++i) {
    comp_index[components[i]] = i;
  }

  ValidateInput();
  
  // Create matrices
  A_burn = BuildMatrix(feed_rate);
  A_off  = BuildMatrix(0.0);
  // Also create exponentiated matrices, assuming constant timestep.
  // Saves repeated cost during timesteps
  EXPA_burn = (A_burn * context()->dt()).exp();
  EXPA_off = (A_off * context()->dt()).exp();

  // Set inventory sizes
  if (compute_startup) EstimateStartup();

  fuel_startup_policy
      .Init(this, &tritium_storage, std::string("Tritium Storage"),
            &fuel_tracker, std::string("ss"),
            startup_inventory, startup_inventory)
      .Set(fuel_incommod, tritium_comp)
      .Start();

  // Tritium Buy Policy Selection:
  if (refuel_mode == "schedule") {
    IntDistribution::Ptr active_dist = FixedIntDist::Ptr(new FixedIntDist(1));
    IntDistribution::Ptr dormant_dist =
        FixedIntDist::Ptr(new FixedIntDist(buy_frequency - 1));
    DoubleDistribution::Ptr size_dist =
        FixedDoubleDist::Ptr(new FixedDoubleDist(1));

    fuel_refill_policy
        .Init(this, &tritium_storage, std::string("Input"), &fuel_tracker,
              buy_quantity, active_dist, dormant_dist, size_dist)
        .Set(fuel_incommod, tritium_comp);

  } else if (refuel_mode == "fill") {

    fuel_refill_policy
        .Init(this, &tritium_storage, std::string("Input"), &fuel_tracker,
              std::string("ss"), reserve_inventory, reserve_inventory)
        .Set(fuel_incommod, tritium_comp);

  } else {
    throw KeyError("Refuel mode " + refuel_mode +
                   " not recognized! Try 'schedule' or 'fill'.");
  }

  tritium_sell_policy.Init(this, &tritium_excess, std::string("Excess Tritium"))
      .Set(fuel_outcommod)
      .Start();
  
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Constructs matrices used for evolving tritium. Input is the rate at which
// tritium is removed from the store and fed into the plasma. This should be
// zero if the plant is switched off. 
Eigen::MatrixXd FlexibleFusionPlant::BuildMatrix(double tritium_consumption_rate) {

  // Add +1 for the plasma
  int N = components.size() + 1;
  
  int plasma = N - 1;
  int storage = comp_index["storage"];
  
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N,N);
      
  // Fill transfer terms
  for (int flow = 0; flow < transfer_rate.size(); flow++) {

    int from = comp_index[transfer_from[flow]];
    int to   = comp_index[transfer_to[flow]];

    double rate = transfer_rate[flow];

    // Off-diagonal gain
    A(to, from) += rate;

    // Diagonal loss to other component
    A(from, from) -= rate;

  }

  // Fill plasma escape terms
  for (int flow = 0; flow < escape_fraction.size(); flow++) {

    int to = comp_index[escape_to[flow]];
      
    double fraction = escape_fraction[flow];

    A(to, plasma) = fraction * (1 - TBE) * tritium_consumption_rate;
	
  }

  // Add constant removal from storage into plasma
  // The indexing is due to the plasma being an inhomogeneous source
  // which does not itself gain tritium. Because the plasma vector element
  // is constant, in the matrix equation, this line corresponds to a constant 
  // removal rate from the storage.
  A(storage, plasma) = -tritium_consumption_rate;
  
  // Add the tritium source term
  // In analogy with the above, this corresponds to a constant production rate
  // in the breeder.
  int breeder = comp_index["breeder"];
  A(breeder, plasma) = TBR * TBE * tritium_consumption_rate;

  // Also add diagonal tritium decay term except in the plasma
  for (int component = 0; component < components.size(); component++) {
    A(component, component) -= pyne::decay_const(tritium_id);
  }

  return A;

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void FlexibleFusionPlant::ValidateInput() {
  
  // Ensure that components contains storage and breeder
  std::vector<std::string> required = {
      "breeder",
      "storage"
  };

  for (auto const& name : required) {
    if (comp_index.count(name) == 0) {
      throw cyclus::ValueError(
          "Required tritium component '" +
          name +
          "' was not defined in input.");
    }
  }

  // Ensure that components DO NOT contain plasma
  if (comp_index.count("plasma") > 0) {
    throw cyclus::ValueError("plasma is a reserved component name");
  }

  // Ensure transfer vectors have the same length
  if (transfer_from.size() != transfer_to.size() ||
      transfer_from.size() != transfer_rate.size()) {
    throw cyclus::ValueError(
        "Transfer vectors must have equal length.");
  }
  
  // Check that transfers are going to/from real places
  for (int flow = 0; flow < transfer_rate.size(); flow++) {

    _require_string(comp_index, transfer_from[flow],
	"Unknown transfer source component: " + 
	transfer_from[flow]);

    _require_string(comp_index, transfer_to[flow],
	"Unknown transfer destination component: " + 
	transfer_to[flow]);

  }
  
  // And escape fractions
  if (escape_to.size() != escape_fraction.size()) {
    throw cyclus::ValueError(
        "Escape vectors must have equal length.");
  }

  // Ensure tritium escape fraction sums to less than one,
  // are all positive, and to real locations
  if (escape_fraction.size() > 0) {
    
    // Check total
    double total = std::accumulate(escape_fraction.begin(),
		    escape_fraction.end(), 0.0);
    if (total > 1.0) {
      throw cyclus::ValueError(
        "Escape fractions must sum to <= 1.");
    }

    // Check positivity
    if (std::any_of(escape_fraction.begin(), escape_fraction.end(),
			    [](double f) {return f < 0;})) {
      throw cyclus::ValueError(
          "All escape fractions must be non-negative.");
    }
    
    // Check validity of destinations
    for (int flow = 0; flow < escape_fraction.size(); flow++) {

      _require_string(comp_index, escape_to[flow],
		      "Unknown escape destination component: " + 
		      escape_to[flow]);
    }

  }
  
  // Ensure startup inventory is greater than reserve  
  if (startup_inventory < reserve_inventory &
		  !compute_startup) {
    throw cyclus::ValueError(
        "Startup inventory must exceed or equal reserve inventory."
	);
  }

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Estimates the startup inventory required by the plant.
// Does so by finding the equilibrium solution for the governing ODE system.
// This amount to solving the linear system when dT/dt = 0
// This cannot estimate reserve inventory as that is determined by tritium loss
// rates from component failures.
// As the storage inventory is fixed, we solve a smaller system of equations
// describing the movement of tritium to all other components.
// We also remove components which are essentially pure sinks, as these do
// not have a true steady-state, other than that set by tritium decay.
void FlexibleFusionPlant::EstimateStartup() {

  double availability = 1.0 - failure_probability;

  Eigen::MatrixXd A_burn_eq = BuildMatrix(feed_rate);
  Eigen::MatrixXd A_off_eq  = BuildMatrix(0.0);
  
  // Setup the linear system, explicitly extracting the source component
  int N = components.size() + 1;
  Eigen::MatrixXd A = availability * A_burn_eq + failure_probability * A_off_eq;
  Eigen::VectorXd Q = A.col(N-1).head(N-1);
  
  int storage = comp_index["storage"];
  // Identify components with no outflow of their own (aside from decay).
  // These are pure sinks -- nothing ever leaves them, so they never reach
  // equilibrium on any realistic timescale, and pre-charging them isn't
  // needed for the plant to sustain operation. Exclude them from the
  // reduced system, just like storage.
  std::vector<double> outflow_rate(components.size(), 0.0);
  for (int flow = 0; flow < transfer_rate.size(); flow++) {
    outflow_rate[comp_index[transfer_from[flow]]] += transfer_rate[flow];
  }

  const double rate_eps = 1e-12;  // effectively-zero transfer rate, in 1/s
  std::vector<bool> excluded(components.size(), false);
  excluded[storage] = true;
  int n_sinks = 0;
  for (int i = 0; i < components.size(); i++) {
    if (i != storage && outflow_rate[i] < rate_eps) {
      excluded[i] = true;
      n_sinks++;
    }
  }

  // Make a reduced matrix since we don't need the storage equation.
  // The storage is given by reserve_inventory. Hence we can solve the smaller
  // system for other components.
  // We also don't need the excluded, pure sink, equations
  int N_red = components.size() - 1 - n_sinks;
  
  // Catch the case where nothing is moving
  if (N_red == 0) {
    startup_inventory = 0.0;
    return;
  }

  Eigen::MatrixXd A_red(N_red, N_red);
  Eigen::VectorXd Q_red(N_red);

  // Populate the matrix, neglecting the unnecessary equations
  int row_idx = 0;
  for (int i = 0; i < N - 1; ++i) {
    if (excluded[i]) continue;

    Q_red(row_idx) = Q(i) + A(i, storage) * reserve_inventory;

    int col_idx = 0;
    for (int j = 0; j < N - 1; ++j) {
      
      if (excluded[j]) continue;

      A_red(row_idx, col_idx) = A(i, j);

      col_idx++;
    }
    row_idx++;
  }

  Eigen::VectorXd x_eq = -A_red.colPivHouseholderQr().solve(Q_red);

  startup_inventory = x_eq.sum() + reserve_inventory;

  // Flag if something has gone wrong
  if ((x_eq.array() < 0.0).any()) {
    throw cyclus::ValueError("Estimated negative tritium densities at equilibrium.\n"
		    "Check transfer/escape values.");
  }

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void FlexibleFusionPlant::Tick() {

  // Use this to check whether tritium is rising or falling
  double previous = tritium_storage.quantity();

  if (ReadyToOperate()) {
    has_started = true;
    fuel_startup_policy.Stop();
    fuel_refill_policy.Start();

    OperateReactor();

  } else {

    // Tritium in the system leaks and decays, but no burning plasma
    OperateReactor(false);

  }

  // Check whether the plant hit or surpassed the reserve inventory
  // Surpassing will have occurred if the current stored tritium is greater
  // than the previous stored tritium after a step
  double current = tritium_storage.quantity();
  bool finished_coasting_down = current <= reserve_inventory || current > previous;
  if (has_started && finished_coasting_down) time_to_sell = true;
  
  // Otherwise, decide whether to send stored tritium to the excess.
  // Requires deciding whether the plant has passed the point at which it
  // should have reached its reserve inventory (if self-sustaining).
  if (time_to_sell && tritium_storage.quantity() > reserve_inventory) {

    // Move tritium in excess of reserve_inventory to the excess
    double transfer_mass = std::max(tritium_storage.quantity() - reserve_inventory, 0.0);
    cyclus::Material::Ptr mat = tritium_storage.Pop(transfer_mass);
    tritium_excess.Push(mat);
  }

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void FlexibleFusionPlant::Tock() {
  RecordInventories(tritium_storage.quantity(), tritium_excess.quantity(),
                    SequesteredTritium());
}

void FlexibleFusionPlant::RecordInventories(double tritium_storage,
                                         double tritium_excess,
                                         double sequestered_tritium) {
  
  auto datum = context()->NewDatum("FFPInventories");
  datum
      ->AddVal("AgentId", id())
      ->AddVal("Time", context()->time())
      ->AddVal("TritiumStorage", tritium_storage)
      ->AddVal("TritiumExcess", tritium_excess)
      ->AddVal("TritiumSequestered", sequestered_tritium);

  // This is necessary because AddVal seems to have a problem
  // being fed a dynamic string.
  column_names_buffer.clear();
  column_names_buffer.reserve(components.size());

  for (int idx = 0; idx < components.size(); idx++) {
    std::string comp_name = components[idx];
    double mass = 0.0;

    if (comp_name == "storage") {
      continue;
    } else if (!tritium_elsewhere[idx].empty()) {
      cyclus::toolkit::MatQuery mq(tritium_elsewhere[idx].Peek());
      mass = mq.mass(tritium_id);
    }

    column_names_buffer.push_back("Tritium" + comp_name);
    datum->AddVal(column_names_buffer.back().c_str(), mass);
  }

  datum->Record();

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double FlexibleFusionPlant::SequesteredTritium() {
  double current_sequestered_tritium = 0.0;

  for (int i = 0; i < components.size(); i++) {
    
    if (tritium_elsewhere[i].empty()
        || i == comp_index["storage"]) {continue;}
    
    cyclus::toolkit::MatQuery mq(tritium_elsewhere[i].Peek());
    current_sequestered_tritium += mq.mass(tritium_id);
  }

  return current_sequestered_tritium;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bool FlexibleFusionPlant::ReadyToOperate() {
  
  // If the plant has failed, increment the recovery counter
  if (has_failed) {

    // Plant restarts
    if (++recovery_counter >= shutdown_duration) {

      recovery_counter = 0;
      has_failed = false;

    // Plant remains off
    } else {
      return false;
    }
  }
	  
  // Determine tritium inventory required to operate
  const double tritium = tritium_storage.quantity();
  if (tritium < startup_inventory && !has_started) {
    return false;
  }
  
  // Check if there is a disruption that prevents operation
  double xi = context()->random_01();
  if (xi < failure_probability) {
    has_failed = true;
    recovery_counter = 0;
    return false;
  }

  // Check if there is enough tritium required to operate: basically operate
  // the reactor and then check if the storage inventory went negative!
  int storage = comp_index["storage"];
  Eigen::VectorXd tritium_vector = EXPA_burn * CurrentTritiumVector();
  if (tritium_vector(storage) < -tritium_eps) return false;

  return true;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// Builds the tritium vector for use with a matrix
Eigen::VectorXd FlexibleFusionPlant::CurrentTritiumVector() {
  int N = components.size() + 1;
  int plasma = N - 1;
  Eigen::VectorXd v(N);
  for (int i = 0; i < N; i++) {
    // Skip the plasma - does not explicitly contain tritium
    // Essentially assumes tritium has zero residence time in
    // the plasma.
    // Contains '1' to act as an inhomogeneous source
    if (i == plasma) {
      v(i) = 1;
    } else if (i == comp_index["storage"]) {
      v(i) = tritium_storage.quantity();
    } else {
      v(i) = tritium_elsewhere[i].quantity();
    }
  }
  return v;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void FlexibleFusionPlant::OperateReactor(bool burn_tritium) {

  double dt = context()->dt();
  
  // Construct tritium vector and evolve it according to 
  // burn rate and transition rates.
  // The final element is the plasma
  Eigen::VectorXd tritium_vector = CurrentTritiumVector();
  std::vector<cyclus::Material::Ptr> popped_mats(components.size());  

  // For a mass balance check, compute the initial mass of tritium.
  // Exclude the '1' plasma term
  double initial_tritium = tritium_vector.sum() - 1.0;

  // Estimate the global tritium post-operation
  double lambda = pyne::decay_const(tritium_id);
  double final_tritium = initial_tritium * std::exp(-lambda * dt);
  if (burn_tritium) {
    double loss = 1.0 - 
  	    std::accumulate(escape_fraction.begin(), escape_fraction.end(), 0.0);
    double constant_gain = burn_rate * (TBR - 1)
	    - feed_rate * loss * (1 - TBE);
    final_tritium += constant_gain / lambda *
	    (1.0 - std::exp(-dt * lambda));
  }
  
  // Choose the matrix based on whether plasma is burning
  // Evolve the densities of tritium in each component
  Eigen::VectorXd new_tritium;
  if (burn_tritium) {
    new_tritium = EXPA_burn * tritium_vector;
  } else {
    new_tritium = EXPA_off * tritium_vector;
  }

  // Update the densities in the appropriate buffers
  double total_tritium = 0.0;
  for (int comp = 0; comp < components.size(); comp++) {
   
    double current_mass = tritium_vector(comp);
    if (new_tritium(comp) < 0.0) {
      cyclus::Warn<cyclus::VALUE_WARNING>("Negative tritium has been produced "
		    "in component " + components[comp] + ": " + 
		    std::to_string(new_tritium(comp)) +"\n"
		    "This will be clamped to zero.");
      new_tritium(comp) = 0.0;
    }
    double new_mass = new_tritium(comp);
    double delta = new_mass - current_mass;

    total_tritium += new_mass;

    // Handle Storage Buffer
    if (comp == comp_index["storage"]) {
      if (delta < -tritium_eps) {
        tritium_storage.Pop(-delta);
      } else if (delta > tritium_eps) {
        tritium_storage.Push(cyclus::Material::Create(this, delta, tritium_comp));
      }
    }
    // Handle other components
    else {
      if (delta < -tritium_eps) {
        tritium_elsewhere[comp].Pop(-delta);
      } else if (delta > tritium_eps) {
        tritium_elsewhere[comp].Push(cyclus::Material::Create(this, delta, tritium_comp));
      }
    }

  }

  // Compare the new masses computed from the matrix solver and the global balance
  if (std::abs(total_tritium - final_tritium) > 1.0e-6) {
      std::cout << "Error in tritium balance: "
                << "matrix solver = " << total_tritium
                << ", global balance = " << final_tritium
                << ", difference = " << (total_tritium - final_tritium)
                << std::endl;
  }

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void FlexibleFusionPlant::_require_string(std::map<std::string, int> string_map,
		std::string required, std::string error_string) {
      
  if (string_map.count(required) == 0) {
    throw cyclus::ValueError(error_string);
  }

}

// WARNING! Do not change the following this function!!! This enables your
// archetype to be dynamically loaded and any alterations will cause your
// archetype to fail.
extern "C" cyclus::Agent* ConstructFlexibleFusionPlant(cyclus::Context* ctx) {
  return new FlexibleFusionPlant(ctx);
}

}  // namespace tricycle
