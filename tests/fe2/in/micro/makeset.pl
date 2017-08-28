#!/usr/bin/perl

$filebase = "pack";
$ext = "st";
$t3d_options = "-d 10. -p 8"; #'-\$ -X';
$decomposition_options = '-sp';
$convert_options = '-v -pr 1.0 1.0 1.0';

$t3d_path = 't3d';
$split_path = 't3d2psifel';
$convert_path = 'con3d++';

sub do_error {
  my $string = shift;
  
  print ($string);
  exit(0);
}

# check for files
if ( !(-e "$filebase"."_0.in.$ext") ) { &do_error("File $filebase"."_0.in.$ext does not exist.\n");}
if ( !(-e "$filebase.out.header") ) { &do_error("File $filebase.out.header does not exist.\n");}
if ( !(-e "$filebase.out.periodic") && ($convert_options =~ /pr/) ) { &do_error("File $filebase.out.periodic does not exist.\n");}

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
      system ("cd $dir; cp -sf ../partitions.$nproc/$filebase.out.* .");
    } else {
      system ("cp $filebase.out $dir/$filebase.out.0");
    }

    print ("Copying header file.\n");
    system ("cp $filebase.out.header $dir/$filebase.out.header");
    system ("cp $filebase.out.periodic $dir/$filebase.out.periodic");
    system ("cd $dir; cp -sf ../normal.in .");

    print ("Converting set...\n");
    system ("$convert_path -f $dir/$filebase -np $nproc $convert_options");
    print ("Finished converting set.\n\n");

    print ("Copying header files...\n");

    for ($i=0; $i<$nproc; $i++) {
	system("cd $dir; \
 cp -sf ../$filebase"."_0.in.$ext $filebase" . "_$i.in.$ext;");
    }

    system("cd $dir; cp -sf ../traction.in ."); 
    print ("\nFinished copying header files.\n\n");
    
    print ("Cleaning up directory...\n");
    system ("rm -rf $dir/*.out.*");
    print ("Finished cleaning up.\n\n");
    
    print ("Finished creating data set for $nproc processors.\n\n");

      
    
  } # end while $nproc = shift
} #end $flag eq '-np'


