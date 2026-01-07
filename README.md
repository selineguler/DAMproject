# DAMproject

present the E-V graphs of MLIPS
present P-V graphs of MLIPS
presenr a-V and c-V of MLIPs (for COF it was relevant, tell if it also relevant of MgO)
<br>

mention your motivation by explaining the Matbench table. There was a lot of data already ,
okay mention about 
<br>

About equation of state studies in elemantarly, we have two types of equation of states which envoles fix symmetry and relaxation EoS.
<br>
With fix symmetry EoS, we fix the symmetry and internal atomic positions of the molecule.
And we compute a single-point energy each time. Applicable structures dependent on the topology, and with the help of symmetry as it already fixes atomic positions.
Relaxed EoS on the otherhand, for each volume, it only fixes the total volume and relax atomic positions and sometimes shape, and computes the minimized energy. This is the correct physical E(V) for materials with internal degrees of freedom.
The point is during the relaxation, atoms move and minimized energy is detected as force goes to 0.
In complex crystals we do relaxation.

I worked on cubic MgO which is a very applicable structure to have.

EXPLAIN WHY RELAXATION IS NOT IMPORTANT WITH MGO. (structure 1/2 1/2 1/2 cubic symmetric etc.)

So MgO was a safe structure to do fix symmetry EoS study to explore MLIPs and also it is well known in literature and won't cause controversy with MLIPs, like it doesn't threat the score of the MLIP, all MLIP should work fine with MgO otherwise it wouldn't make it to the benchmark anyways.




In the beginning I wanted to actually somehow reproduce the Matbench discovery parameters, not exactly reproduce but at least focus on how they are concluded, but these are all evaluated over tens or hundreds of thousands of structures and via studying this simple structure, my results mostly focuses on simply relating to these parameters. It is a likelihood-based mini-mini benchmark, cause I neither used the whole MLIPs nor enough structures to make a reliable analysis, but the project in general lead me to understand more on MLIPs and these machine-learning testing parameters that are important to choose of a calculator for a possible future system to work on.

Metrics I focused on for this MgO-centric benchmark is MAE and RMSE ***insert formulas here***


I also implemented log-likelihood to score the MLIPs.

We do the fitting with EOS but EOS is literally fundemental physics so likelihood must be scoring the MLIP (right?)

We have reference energies, and calculated energies; and we assume MLIP errors is a function with some uncertainty and that function is mostly Gaussian.
Then we do the likelihood with the prev uncertainty and ref and calculated energies.
And conclude that larger log-likelihood => smaller squared error => better MLIP for the MgO.
You can do the same for pressures and explain why you can do the same for pressures.
At the end how we will be able to relate this to the matbench leaderboard?
We cannot, that failed, matbench metrics are way too global and estimated with way more structures so for a single MgO we cannot say the likelihood we calculated for CPS corresponds to the one on the benchmark. There is not even max we can do, we literally cannot associate our likelihood calculations to the benchmark anyhow. Only compare the calculators with the literature and if the order is same within the discovery benchmark it's good, but even if it's not it doesn't indicate anything is wrong because you literally worked with a single molecule. But at least you know the core idea of how these MLIPs are benchmarked. That is the motivation of the project.




