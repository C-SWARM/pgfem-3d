#!/usr/bin/perl

use warnings;
if ($#ARGV==-1){
        print "usage: ./local_makeset.pl -np [NP] [el_size] [el_order]\n";
        print "\t[NP]       : number of processors\n";        
	print "\t[el_size]  : element size\n";
        print "\t[el_order] : element order linear 1, quadratic 2\n";
	exit 0;
}

$filebase = "beam";
$ext = "st";
$flag = shift;

if ( $flag eq 'clean') {
    system ("rm -rf $filebase.out *CPU partitions.*");
    exit 1;
}

# read inputs
$nproc = $ARGV[0];
$el_size = $ARGV[1];
$el_order = $ARGV[2];

# t3d/ t3d2psifel/ con3d++ options
$t3d_options = "-d $el_size -k $el_order -r 1 -p 8"; #'-\$ -X';
$decomposition_options = '-sp';
$convert_options = '-v -HumanReadable'; #'-pr 200.0 200.0 200.0';

# path settings
$t3d_path = '/afs/crc.nd.edu/group/cswarm/bin/t3d';
$split_path = '/afs/crc.nd.edu/group/cswarm/bin/t3d2psifel';
$convert_path ='/afs/crc.nd.edu/group/cswarm/bin/con3d++';

sub do_error {
my $string = shift;
 
print ($string);
exit(0);
}

# check for files
if ( !(-e "$filebase.out.periodic") && ($convert_options =~ /pr/) ) { &do_error("File $filebase.out.periodic does not exist.\n");}

    print ("Preparing data set for $nproc processors.\n\n");

    if ( !(-e "$filebase.out") ) {
      print "Meshing data set.\n\n";
      
      system "$t3d_path -i $filebase.t3d -o $filebase.out $t3d_options";
      system "cp $filebase.out $filebase.out ";
    }
    
    $dir = "$filebase"."_$nproc" . "CPU";

    mkdir ($dir);

    if ( $nproc != 1 && !(-e "partitions.$nproc") ) {
      print ("Splitting file...\n");
      mkdir ("partitions.$nproc");
      system ("$split_path -np $nproc $filebase.out partitions.$nproc/$filebase.out $decomposition_options");
      print ("Finished splitting file.\n\n");
    } 

    if ( $nproc != 1 ) {
      system ("cd $dir; cp -sf ../partitions.$nproc/$filebase.out.* .");
    } else {
      system ("cp $filebase.out $dir/$filebase.out.0");
    }

    print ("Converting set...\n");
    system ("$convert_path -f $dir/$filebase -np $nproc $convert_options -d $dir");

    print ("Finished converting set.\n\n");

    print ("Copying in model_params.in\n\n");
    system ("cp model_params.in $dir");
    print ("Copying header files...\n");
    for ($i=1; $i<$nproc; $i++) {
	system("cd $dir; \
 cp -sf $filebase"."_0.in.$ext $filebase" . "_$i.in.$ext" );
    }
    system("cd $dir; cp -sf ../traction.in ."); 
    print ("\nFinished copying header files.\n\n");
    print ("Cleaning up directory...\n");
    system ("rm -rf $dir/*.out.*");
    print ("Finished cleaning up.\n\n");
    print ("Finished creating data set for $nproc processors.\n\n");
