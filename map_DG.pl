#!/usr/bin/perl
# 
#   *************
#   * map_DG.pl *
#   *************
#
# This soft compute a rough estimate of a free energy landscape using
# two arbitrary structural coordinates. Both of them can came, for
# example, from a molecular dynamics trajectory. This computation is 
# made through Boltzmann dependance between population and Gibbs free
# energy. We aware the reader that, without reaching the system 
# ergodicity, this computation is only an estimate and should be seen
# more qualitatively than quantitatively.
#
# Usage:
# % map_DG.pl file
# 
# Where file is a text file which contains 3 columns. These are, 
# respectively, an index number, the X and Y values. To get this file,
# it might be useful for the use the 'extract_column.pl' software which
# selects the desired columns from a more extended file and create a new
# one. 
# The software will first identified the minimum and maximum values on
# each dimensions. The user will then have to choose between using these
# values or defining new ones.The user will be prompted to give the number
# of wanted bins. This value is identical for both dimension but may be
# adapted if someone ask. For the image generation, made with gnuplot, 
# the user will have to provide labels for the reaction coordinates.
#
# The software will produced three files:
# + map.dat
# + map.gnu
# + map.png
# 'map.dat' contains all the data whereas 'map.gnu' is a script used
# by gnuplot for generating the landscape image such as 'map.png'.
# 
# Enjoy, Florent Barbault.
#
$|;
use strict;
use Math::Trig;
#---------------------------------------------------------------------
# Variable declaration
#
my $file=$ARGV[0];
my $i=0;
my $j=0;
my @tab=();
my @table=();
my @X=();
my @Y=();
my @xx=();
my @yy=();
my @oc=();
my $tot=0;
my $xmin=0;
my $ymin=0;
my $xmax=0;
my $ymax=0;
my $bin_tot=0;
my $dx=0;
my $dy=0;
my $R=8.314;
my $T=300;                #Temperature, this may be changed if needed
my $kcal=4180;
my $factor=0;
my $occur=0;
my $xcoord=0;
my $ycoord=0;
my $xbas=0;
my $xhaut=0;
my $ybas=0;
my $yhaut=0;
my $DG=0;
my $total=0;
my $OCCUR=0;
my $OCCUR_max=0;
my $max=0;
my $ind=0;
my $choix="";
my $xlabel="";
my $ylabel="";
my $bin_grid=0;
my $ask_gnu="";
#---------------------------------------------------------------------
# Reading data
#
open (FILE, "< $file");
while (<FILE>)
{
  @tab = split(' ',$_);
  $X[$i]=$tab[1];
  $Y[$i]=$tab[2];
  $i=$i+1;
}
$tot=$i;
print "=> You have $tot data\n";
#---------------------------------------------------------------------
# Getting the lowest and highest values of x and y
#
$xmin=$X[0];
$xmax=$X[0];
$ymin=$Y[0];
$ymax=$Y[0];
for ($i=0;$i<$tot;$i++)
{
 if ($X[$i] < $xmin){ $xmin=$X[$i];}
 if ($X[$i] > $xmax){ $xmax=$X[$i];}
 if ($Y[$i] < $ymin){ $ymin=$Y[$i];}
 if ($Y[$i] > $ymax){ $ymax=$Y[$i];}
}
print "=> your values are ranked between\n";
print "      $xmin < X < $xmax\n";
print "      $ymin < Y < $ymax\n";
#---------------------------------------------------------------------
# Ask if the range is ok or provide other values
#
print "\nDo you want to keep these extremum values as axis limit?";
print "\nType Y for yes, all other values will lead to no\t";
$choix=<STDIN>;
chop($choix);
if ($choix eq 'Y')
 {

 } 	
else 
 {
   print "\nX min value:\t";
   $xmin=<STDIN>;
   chop($xmin);
   print "X max value:\t";
   $xmax=<STDIN>;
   chop($xmax);
   print "Y min value:\t";
   $ymin=<STDIN>;
   chop($ymin);
   print "Y max value:\t";
   $ymax=<STDIN>;
   chop($ymax);
}
#---------------------------------------------------------------------
# Ask how many bins the user want for X and Y
#
print "=> How many bin you want in the X and Y directions:  ";
$bin_tot=<STDIN>;
chop($bin_tot);
#---------------------------------------------------------------------
# Ask for X and Y labels
#
print "=> Label for the reaction coordinate in the X directions:  ";
$xlabel=<STDIN>;
chop($xlabel);
print "=> Label for the reaction coordinate in the Y directions:  ";
$ylabel=<STDIN>;
chop($ylabel);
#---------------------------------------------------------------------
# Delineates the areas and count the occurence
#
$dx=($xmax-$xmin)/$bin_tot;
$dy=($ymax-$ymin)/$bin_tot;
open (FOUT, "> temp.dat");
open (GOUT, "> map.gnu");
for ($i=1;$i<=$bin_tot;$i++)
{
 for ($j=1;$j<=$bin_tot;$j++)
 {
   $xcoord=$xmin+$i*$dx-$dx/2;
   $ycoord=$ymin+$j*$dy-$dy/2;
   $occur=0;
   for ($ind=0;$ind<$tot;$ind++)
   {
     $xbas=$xcoord - $dx/2;
     $xhaut=$xcoord + $dx/2;
     $ybas=$ycoord - $dy/2;
     $yhaut=$ycoord + $dy/2;
     if ($X[$ind] >= $xbas)
      {
       if ($X[$ind] < $xhaut)
        {
         if ($Y[$ind] >= $ybas)
          {
           if ($Y[$ind] < $yhaut)
            {
              $occur=$occur+1;
            }
          }
        }
      }
   }
   printf FOUT ("\n %6.2f   %6.2f   %6.2f",$xcoord,$ycoord,$occur);
 }
}
close(FOUT);
#---------------------------------------------------------------------
# Transform occurences in DG and make the map.dat file
#
$factor=$R*$T/$kcal;
open (FG, "> map.dat");
$i=0;
$total=0;
$OCCUR=0;
$OCCUR_max=0;
$max=0;
open (FILE, "< temp.dat");
while (<FILE>)
{
  @table = split(' ',$_);
  $xx[$i]=$table[0];
  $yy[$i]=$table[1];
  $oc[$i]=$table[2];
  $OCCUR=$OCCUR+$oc[$i];
  if ($oc[$i] > $OCCUR_max){$OCCUR_max=$oc[$i];}
  $i=$i+1;
}
$total=$i;
for ($ind=0;$ind<$total;$ind++)
{
 if ($oc[$ind] != 0.0){ $DG=sqrt((-$factor*log($oc[$ind]/$OCCUR_max)*(-$factor*log($oc[$ind]/$OCCUR_max)))); if ($DG > $max){$max=$DG;}}
 if ($oc[$ind] == 0.0){ $DG=0;}
 printf FG ("\n %6.2f   %6.2f   %6.2f",$xx[$ind],$yy[$ind],$DG);
}
close(FG);
close(FILE);
system ("rm temp.dat");
$bin_grid=2*$bin_tot;
#---------------------------------------------------------------------
# Create GNUPLOT file 'map.gnu' 
#
print GOUT "set xrange [$xmin:$xmax]\n";
print GOUT "set yrange [$ymin:$ymax]\n";
print GOUT "set dgrid3d $bin_tot,$bin_tot splines\n";
print GOUT "set xlabel '$xlabel'\n";
print GOUT "set ylabel '$ylabel'\n";
print GOUT "set pm3d at b\n";
print GOUT "unset surface\n";
print GOUT "unset key\n";
print GOUT "unset clabel\n";
print GOUT "set size square\n";
print GOUT "set view map\n";
print GOUT "set palette rgbformulae -33,-13,-10\n";
print GOUT "set term png enhanced size 1000,1000 font \"arial,20\"\n";
print GOUT "set output \"map.png\"\n";
print GOUT "splot \"map.dat\" w li\n";
print GOUT "quit\n";
close (GOUT);
#---------------------------------------------------------------------
# Ask for the generation of gnuplot to get 'map.png'
#
print "\nDo you want to generate a map image with gnuplot?";
print "\nType Y for yes, all other values will lead to no\t";
$ask_gnu=<STDIN>;
chop($ask_gnu);
if ($ask_gnu eq 'Y')
 {
   system ("gnuplot < map.gnu");
 } 	
else 
 {
   
 }
#---------------------------------------------------------------------
# End
#
print "\n The end, enjoy\n";


