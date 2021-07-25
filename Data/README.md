Data used for this vaccine allocation study can be downloaded at https://bit.ly/3zexGiu

File description:

 * __data_node.txt__ contains information on nodes of the rasterized (50 km x 50 km) world network. Each row corresponds to a node and contains the following information: node ID, country (ISO-2 code), population, followed by a list of node IDs that represent nodes in the neighborhood of this node.

 * __data_dm.txt__ contains information on decision-making agents which represents countries. Each row corresponds to a country and contains the following information: country ID (ISO-2 code), country population, initial susceptible population, initial infectious population, disease transmission rate (beta), disease recovery rate (gamma), per-period vaccination capacity, followed by a list of node IDs that represent nodes controlled by this country.

 * __data_mobility.txt__ contains information on the human mobility flows in the rasterized world network. Each row is a link (i,j) in the world network and contains the following information: node ID of head node (i), node ID of tail node (j), mobility rate (p_ij) and a binary value indicating whether nodes i and j are controlled by the same country (1) or not (0).
