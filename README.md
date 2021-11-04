# map_DG
calculate free energy landscape from 2 metrics

This soft compute a rough estimate of a free energy landscape using two arbitrary structural coordinates. Both of them can came, for example, from a molecular dynamics trajectory. This computation is made through Boltzmann dependance between population and Gibbs free energy. We aware the reader that, without reaching the system ergodicity, this computation is only an estimate and should be seen more qualitatively than quantitatively.

Usage:
% map_DG.pl file

Where file is a text file which contains 3 columns. These are, respectively, an index number, the X and Y values. To get this file, it might be useful for the use the 'extract_column.pl' software which selects the desired columns from a more extended file and create a new one. The software will first identified the minimum and maximum values on each dimensions. The user will then have to choose between using these values or defining new ones.The user will be prompted to give the number of wanted bins. This value is identical for both dimension but may be adapted if someone ask. For the image generation, made with gnuplot, the user will have to provide labels for the reaction coordinates.

The software will produced three files:
 + map.dat
 + map.gnu
 + map.png

'map.dat' contains all the data whereas 'map.gnu' is a script used by gnuplot for generating the landscape image such as 'map.png'.

Enjoy, Florent Barbault.
If this software is useful for you, please cite:
"Molecular Dynamics Simulation of a RNA Aptasensor" M. Ruan, M. Seydou, V. Noël, B. Piro, F. Maurel and F. Barbault*, Journal of Physical Chemistry B, 2017 (16) 4071-80 
