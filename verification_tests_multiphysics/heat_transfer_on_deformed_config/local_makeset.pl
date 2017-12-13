#!/usr/bin/perl

$filebase = "box";
$ext = "st";
$t3d_options = "-d 0.1 -p 8"; #'-\$ -X';
$decomposition_options = '-sp -rn';
$convert_options = '-v -HumanReadable'; 

$t3d_path = '$CSWARM_GROUP/t3d';
$split_path = '$CSWARM_GROUP/t3d2psifel';
$convert_path = 'con3d++';

sub do_error {
  my $string = shift;
  
  print ($string);
  exit(0);
}

# check for files
if ( !(-e "$filebase.out.periodic") && ($convert_options =~ /pr/) ) { &do_error("File $filebase.out.periodic does not exist.\n");}

$flag = shift;

if ( $flag eq '-np' ) {
  # iterate through every number specified
  while ( $nproc = shift ) {
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

#    system("cd $dir; cp -sf ../traction.in ."); 
    print ("\nFinished copying header files.\n\n");
    
    print ("Cleaning up directory...\n");
    system ("rm -rf $dir/*.out.*");
    print ("Finished cleaning up.\n\n");
    
    print ("Finished creating data set for $nproc processors.\n\n");

      
    
  } # end while $nproc = shift
} #end $flag eq '-np'
