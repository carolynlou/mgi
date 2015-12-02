#!/usr/bin/perl

use strict;
use warnings;

use ClusterModule;
use FileHandle;
use IO::File;
use Getopt::Long;


my $object= new ClusterModule();
#my $object = new ClusterModule(~/../ascott/Projects/HotSpot3D/Comparison/ClusterLikeSpacePAC/chosen.cancer.heteromer.hup, ~/../ascott/Projects/HotSpot3D/Comparison/ClusterLikeSpacePAC/Heteromer/NewPValuePairwise/chosen.cancer.heteromer, ~/../ascott/Projects/HotSpot3D/Comparison/ClusterLikeSpacePAC/single_molecule_pvalue.pairwise, ~/../ascott/Projects/HotSpot3D/Comparison/ClusterLikeSpacePAC/Heteromer/NewPValueClusters/chosen.cancer.heteromer, ~/../ascott/Projects/HotSpot3D/Comparison/ClusterLikeSpacePAC/single_molecule_pvalue.intra.20..05.10.clusters, ~/../dinglab/medseq/AA_Conversion/AAabbrv, ~/../dinglab/medseq/Preprocessing_Output_20141023/proximityFiles/cosmicanno/, 1000000, 1>heteromer.pvalues.intra.20..05.10.1M.out, 2>heteromer.pvalues.intra.20..05.10.1M.err);

$object->process();
print $object->getDistances();
print "Hello";
