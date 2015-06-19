#!/usr/bin/perl

$filebase = "box_LT";
$ext = "st";
$t3d_options = "-d .1 -p 8 "; #'-\$ -X';
$decomposition_options = '-sp -rn';
$convert_options = '-v'; #'-pr 200.0 200.0 200.0';

$t3d_path = '/cswarm/tools/bin/t3d';
$split_path = '/cswarm/tools/bin/t3d2psifel';
$convert_path = '/cswarm/tools/bin/con3d++';

sub do_error {
  my $string = shift;
  
  print ($string);
  exit(0);
}

# check for files
if ( !(-e "$filebase"."_0.in.$ext") ) {
    &do_error("File $filebase"."_0.in.$ext does not exist.\n");
}
if ( !(-e "$filebase.out.header") ) {
    &do_error("File $filebase.out.header does not exist.\n");
}
if ( !(-e "$filebase.out.periodic") && ($convert_options =~ /pr/) ) {
    &do_error("File $filebase.out.periodic does not exist.\n");
}

$flag = shift;

if ( $flag eq 'clean') {
  system ("rm -rf $filebase.out *CPU partitions.*");
  $flag = shift;
}

if ( $flag eq '-np' ) {

  # iterate through every number specified
  while ( $nproc = shift ) {
    print ("Preparing data set for $nproc processors.\n\n");

    if ( !(-e "$filebase.out") ) {
      print "Meshing data set.\n\n";
      
      system "$t3d_path -i $filebase.t3d -o $filebase.out $t3d_options";
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
      system ("cd $dir; cp -s ../partitions.$nproc/$filebase.out.* .");
    } else {
      system ("cp $filebase.out $dir/$filebase.out");
    }

    print ("Copying header file.\n");
    system ("cp $filebase.out.header $dir/$filebase.out.header");

    print ("Converting set...\n");
    system ("$convert_path -f $dir/$filebase -np $nproc $convert_options");
    print ("Finished converting set.\n\n");

    print ("Copying header files...\n");
    print ("[" . (' ' x $nproc) . "]");


    for ($i=0; $i<$nproc; $i++) {
      system("cd $dir; cp -s ../$filebase"."_0.in.$ext $filebase" . "_$i.in.$ext");
    }
    print ("\nFinished copying header files.\n\n");
    
    print ("Cleaning up directory...\n");
#    system ("rm -rf partitions.*");
    system ("rm -rf $dir/*.out.*");
    print ("Finished cleaning up.\n\n");
    
    print ("Finished creating data set for $nproc processors.\n\n");

      
    
  } # end while $nproc = shift
} #end $flag eq '-np'


