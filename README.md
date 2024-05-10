## Identifying Peer Effects on Student Academic Effort
This is a replication code of "Identifying Peer Effects on Student Academic Effort" by Aristide Houndetoungan, Cristelle Kouame, and Michael Vlassopoulos.

- Files `SourceR.R` and `SourceCpp.cpp` contain functions built in **R** and **cpp** to be used in other files.
- `0_Inschool.do` extracts the part of the data set to be used from the Add Health data set.
- `0A_data.R` prepares the data set to be used.
- `1_exogenous_network.R` replicates the peer effect model estimation assuming that the network is exogenous
- `2_network_formation.R` estimates the network formation model.
- `3_endogenous_network.R` replicates the peer effect model estimation controlling for network endogeneity.
- `4_shocks.R` simulates shocks on alpha_s and c_s.
- `4_montecarlo.R` replicates our simulation study.
