Data used for this vaccine allocation study can be downloaded at https://bit.ly/3zexGiu

File description:

 * __data_node.txt__ contains information on nodes of the rasterized (50 km x 50 km) world network. Each row corresponds to a node and contains the following information: node ID, country (ISO-2 code), population, followed by a list of node IDs that represent nodes in the neighborhood of this node.

 * __data_dm.txt__ contains information on decision-making agents which represent countries. Each row corresponds to a country and contains the following information: country ID (ISO-2 code), country population (P_i), initial susceptible population (S_i(O)), initial infectious population (I_i(O)), disease transmission rate (beta), disease recovery rate (gamma), per-period vaccination capacity (Gamma_k), followed by a list of node IDs that represent nodes controlled by this country.

 * __data_mobility.txt__ contains information on the human mobility flows in the rasterized world network. Each row is a link (i,j) in the world network and contains the following information: node ID of head node (i), node ID of tail node (j), mobility rate (p_ij) and a binary value indicating whether nodes i and j are controlled by the same country (1) or not (0).

* The folder __scenario_data__ contains 300 randomly generated scenarios used in the study. Scenario files are named in the format __X104_EP_ID.txt__ where EP represents the level of uncertainty in the vaccine administration process (epsilon) written in percent and ID is the identifier of this scenario. Three levels of uncertainty are considered (10%, 20% and 30%) and 100 random scenarios are available for each level of uncertainty. In each scenario file, the top row contains the ID of all nodes of the world network. The second row contains node-based mean vaccine efficiency rates (theta_i) which are unknown to decision-making agents (countries); and the next 104 rows contain node-based observed mean vaccine efficiency rates corresponding to realizations of the random variables (theta_i(t)) for all 104 time periods modeled in this study.

The file Data_countries.pdf contains a summary of country-based data used in this study.
