# Molecular Dynamics Simulation of Lipoteichoic Acid Synthase LtaS from Listeria monocytogenes (PDB ID: 4UOR)

## Summary 
Molecular dynamics (MD) simulation is a powerful computational technique used to study the dynamic behavior of biomolecules at the atomic level. In this report, we performed MD simulations on the lipoteichoic acid synthase LtaS from Listeria monocytogenes in complex with glycerol phosphate, using the protein structure with PDB ID 4UOR. The simulations were conducted on a GPU-enabled Google Colab platform.

## Methods:
To prepare the system for simulation, we employed the following minimization and equilibration steps:

### Energy Minimization:

- We used the steepest descent minimization algorithm to relax the system.
- The minimization process was terminated when the maximum force reached a threshold of 1000 kJ/mol/nm.
- The minimization step size was set to 0.01, and a maximum of 50,000 steps were performed.

### Ionization:

- This step involved adding ions to neutralize the system.
- The same minimization algorithm and parameters as in the previous step were used.

### NVT Equilibration:

- The protein was restrained using position restraints.
- MD simulations were performed using the leap-frog integrator for 50,000 steps (100 ps).
- The time step was set to 0.002 ps (2 fs).
- Energies, logs, and compressed coordinates were saved every 1.0 ps.

## NPT Equilibration:

- The system was equilibrated under constant pressure and temperature conditions.
- Similar simulation parameters as in the NVT equilibration step were used.

### Production MD Simulation:

- The equilibrated system was subjected to a production MD simulation.
- The simulation was run for 50,000 steps (100 ps) with a time step of 0.002 ps (2 fs).
- Energies, logs, and compressed coordinates were saved every 2.0 ps.

## Results:
The MD simulations provided insights into the dynamic behavior of the lipoteichoic acid synthase LtaS in complex with glycerol phosphate. The trajectory data generated during the simulations can be analyzed to study the conformational changes, stability, and interactions of the protein-ligand complex.

### Loaded Cleaned Protein:

![Cleaned Protein](https://github.com/paulshamrat/230427_LtaS/blob/main/results/cleaned_protein.png)

The initial step involves loading the protein structure (4UOR.pdb) into the simulation environment.
Before simulation, it's common to perform cleaning procedures to remove any artifacts or inconsistencies in the protein structure, such as missing atoms or incorrect bond lengths.
By loading the cleaned protein, you ensure that the simulation starts with a reliable and accurate representation of the protein.

### Solvated System:

![Solvated System](https://github.com/paulshamrat/230427_LtaS/blob/main/results/solvated_system.png)

After loading the cleaned protein, the system is solvated by adding water molecules around the protein.
Solvating the system is crucial to simulate the protein's behavior in an aqueous environment, mimicking its natural surroundings.
The number of water molecules and the size of the simulation box are typically determined based on the desired concentration and system stability.

### Potential Energy during Minimization:

![Potential Energy](https://github.com/paulshamrat/230427_LtaS/blob/main/results/potential_energy.png)

Before running molecular dynamics simulations, it's essential to minimize the system's potential energy to remove any steric clashes or high-energy configurations.
During the energy minimization process, the potential energy of the system is iteratively reduced by adjusting the atomic coordinates.
Monitoring the potential energy helps ensure that the system reaches a stable energy state before proceeding to the equilibration phase.

### Temperature during 1000 ps Equilibration (NVT):

![Temperature](https://github.com/paulshamrat/230427_LtaS/blob/main/results/temperature.png)

The equilibration phase (NVT ensemble) aims to bring the system to the desired temperature and allow it to adjust to the thermal environment.
Monitoring the temperature during equilibration ensures that the system reaches and maintains the target temperature.
Deviations from the target temperature can indicate issues with the simulation setup or energy transfer within the system.

### Pressure during 1000 ps Equilibration (NPT):

![Pressure](https://github.com/paulshamrat/230427_LtaS/blob/main/results/pressure.png)

The NPT ensemble includes both temperature and pressure control, allowing the system to adjust to the desired pressure conditions.
Monitoring the pressure during equilibration ensures that the system reaches and maintains the target pressure.
Deviations from the target pressure can indicate issues with the simulation setup, system density, or interactions between the protein and solvent.

### Density during 1000 ps Equilibration (NPT):

![Density](https://github.com/paulshamrat/230427_LtaS/blob/main/results/density.png)

Density refers to the mass per unit volume and provides information about the packing and compactness of the system.
Monitoring the density during equilibration helps ensure that the system reaches a stable density consistent with the experimental conditions or desired simulation setup.
Changes in density may indicate system expansion or contraction, which can impact the protein's interactions and behavior.

### RMSD (Root-Mean-Square Deviation) Plot:

![RMSD Plot](https://github.com/paulshamrat/230427_LtaS/blob/main/results/rmsd_plot.png)

The RMSD plot shows the deviation of the protein's structure from its initial reference structure over time.
It provides insights into the stability and structural changes of the protein during the MD simulation.
Large fluctuations in RMSD may indicate conformational changes or structural rearrangements.

### RMSF (Root-Mean-Square Fluctuation) Plot:

![RMSF CA Plot](https://github.com/paulshamrat/230427_LtaS/blob/main/results/rmsf_ca_plot.png)

The RMSF plot displays the average atomic fluctuations of the protein residues during the MD simulation.
It helps identify flexible or highly mobile regions within the protein.
Higher RMSF values suggest greater flexibility, while lower values indicate more rigid regions.

### Rg (Radius of Gyration) Plot:

![Rg Plot](https://github.com/paulshamrat/230427_LtaS/blob/main/results/rg_plot.png)

The Rg plot shows the average distance of the protein atoms from the center of mass over time.
It provides information about the compactness or expansion of the protein during the simulation.
An increasing Rg value may indicate protein unfolding or structural expansion.

### RMSD Average Linkage Hierarchical Clustering:

![RMSD Average Linkage Hierarchical Clustering](https://github.com/paulshamrat/230427_LtaS/blob/main/results/RMSD_Average_linkage_hierarchical_clustering.png)

Hierarchical clustering is a method to group similar conformations based on their RMSD values.
It helps identify distinct clusters or conformational states of the protein during the simulation.
The clustering dendrogram or heat map can reveal different structural ensembles or transitions.

### Cartesian Coordinate PCA (Principal Component Analysis):

![Cartesian Coordinate PCA](https://github.com/paulshamrat/230427_LtaS/blob/main/results/Cartesian_coordinate_PCA.png)


PCA is a dimensionality reduction technique that captures the dominant motions of the protein.
The Cartesian Coordinate PCA plot displays the projection of the protein's motion onto the principal components.
It helps identify collective motions, such as hinge motions or domain movements, that contribute most to the conformational changes.

### Pairwise Distance PCA:

![Pairwise Distance PCA](https://github.com/paulshamrat/230427_LtaS/blob/main/results/Pairwise_distance_PCA.png)

Pairwise Distance PCA analyzes the distances between pairs of residues during the MD simulation.
It provides insights into the dynamic interactions and rearrangements of specific regions within the protein.
The PCA plot based on pairwise distances helps visualize the correlated motions of residue pairs.

### Total SASA (Solvent-Accessible Surface Area) Plot:

![Total SASA](https://github.com/paulshamrat/230427_LtaS/blob/main/results/Total_SASA.png)

The Total SASA plot displays the time evolution of the protein's exposed surface area.
It can indicate changes in protein solvation, solvent accessibility, or exposure of specific regions.
Fluctuations in SASA can suggest conformational changes or interactions with solvent molecules.

### SASA Autocorrelation Plot:

![SASA Autocorrelation](https://github.com/paulshamrat/230427_LtaS/blob/main/results/sasa_autocorrelation.png)

The SASA autocorrelation plot measures the correlation of the SASA values at different time intervals.
It helps identify the timescale at which the protein's surface area fluctuations persist or decay.
Autocorrelation analysis provides insights into the timescales of protein motions and relaxation.


## Conclusion:
In this case study, we performed MD simulations on the lipoteichoic acid synthase LtaS from Listeria monocytogenes in complex with glycerol phosphate. The simulations provided valuable information about the dynamics of the system and can serve as a basis for further analysis and interpretation of the protein-ligand interactions. The computational approach used in this study demonstrates the utility of MD simulations in understanding the behavior of biomolecules at the atomic level.
