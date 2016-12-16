#!/usr/bin/perl

sub do_error {
  my $string = shift;

  print ($string);
  exit(0);
}

#$t3d_path = '/cswarm/tools/bin/t3d';
#$split_path = '/cswarm/tools/bin/t3d2psifel';
#$convert_path = '/cswarm/tools/bin/con3d++';

$t3d_path = '$CSWARM_GROUP/t3d';
$split_path = '$CSWARM_GROUP/t3d2psifel';
$convert_path = '$CSWARM_GROUP/con3d++';

$var_clean = 0;
$nproc = 0;
$var_d = 1.0;
$var_k = 1;
$var_pr = "";
$filebase = "";
$continue = 1;

while($continue eq 1){
  $flag = shift;

  if ( $flag eq 'clean') {
    $var_clean = 1;
  }

  if( $flag eq '-f'){
    $filebase = shift;
  }

  if( $flag eq '-np'){
    $nproc = shift;
  }

  if( $flag eq '-d'){
    $var_d = shift;
  }

  if( $flag eq '-k'){
    $var_k = shift;
  }

  if( $flag eq '-pr'){
    $temp = shift;
    $var_pr = "-pr $temp";
  }

  if( $flag eq '-h'){
    print("======================================================\n");
    print("usage of makeset.pl\n");
    print("-f filebase_name : filebase \n");
    print("-np # : number of process to decompose \n");
    print("-d # : mesh size \n");
    print("-k # : order of element \n");
    print("-pr # # #: size of preodic cell \n");
    print("clean : clean up \n");
  }
  
  if($flag eq ""){
     $continue = -1;
  }
}

if($var_clean eq 1){
  print("\n");
  print("clean up before proceed : rm -rf *.out *CPU partitions.* \n");
  system ("rm -rf *.out *CPU partitions.*");
}

$ext = "st";
$t3d_options = "-d $var_d -p 8 -k $var_k "; #'-\$ -X';
$decomposition_options = '-nt 12 -sp';
#$decomposition_options = '-sp';
$convert_options = "-v $var_pr";

if( $nproc gt 0)
{
	
  print("\n");
  print("======================================================\n");
  print("makeset.pl options \n");
  print("filebase_name: \t\t\t $filebase \n");
  print("number of process to decompose:  $nproc \n");
  print("mesh size: \t\t\t $var_d \n");
  print("order of element: \t\t $var_k \n");
  print("preodic conditions: \t\t $var_pr \n");

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

  print("start convert ...\n");
  print ("Preparing data set for $nproc processors.\n\n");

  if ( !(-e "$filebase.out") ) {
    print "Meshing data set.\n\n";
    print ("$t3d_path -i $filebase.t3d -o $filebase.out $t3d_options \n");
    system "$t3d_path -i $filebase.t3d -o $filebase.out $t3d_options";
  }
    
  $dir = "$filebase"."_$nproc" . "CPU";
  mkdir ($dir);

  if ( $nproc != 1 && !(-e "partitions.$nproc") ) {
    print ("Splitting file...\n");
    print ("$split_path -np $nproc $filebase.out partitions.$nproc/$filebase.out $decomposition_options \n");
    mkdir ("partitions.$nproc");
    system ("$split_path -np $nproc $filebase.out partitions.$nproc/$filebase.out $decomposition_options");
    print ("Finished splitting file.\n\n");
  }

  if ( $nproc != 1 ) {
    if ( $nproc < 32768) {
      system ("cd $dir; cp -s ../partitions.$nproc/$filebase.out.* .");
    } else {
      for ($ia=0; $ia<10; $ia++)
      {
        system ("cd $dir; cp -s ../partitions.$nproc/$filebase.out.*$ia .");
      }
    }
  } else {
    system ("cp $filebase.out $dir/$filebase.out.0");
  }

  print ("Copying header file.\n");
  system ("cp $filebase.out.header $dir/$filebase.out.header");
  
  if ( ($convert_options =~ /pr/) ) {
    print ("Copying periodic file.\n");
    system ("cp $filebase.out.periodic $dir/$filebase.out.periodic");  
  }

  print ("Converting set...\n");
  system ("$convert_path -f $dir/$filebase -np $nproc $convert_options");
  print ("Finished converting set.\n\n");

  print ("Copying solver files...\n");
  for ($i=0; $i<$nproc; $i++) {
    system("cd $dir; cp -sf ../$filebase"."_0.in.$ext $filebase" . "_$i.in.$ext");
  }

  print ("\nFinished copying solver files.\n\n");
    
  print ("Cleaning up directory...\n");
  if ( $nproc < 32768) {
    system ("rm -rf $dir/*.out.*");
  } else {
    for ($ia=0; $ia<10; $ia++)
    {
      system ("rm -rf $dir/*.out.*$ia");
    }
    system ("rm -rf $dir/*.out.*");
  }
  print ("Finished cleaning up.\n\n");

  # link normal.in file
  if ( -e "normal.in" ) {
    system("cd $dir; cp -sf ../normal.in .");
  }

  # link model_params.in file
  if ( -e "model_params.in" ) {
    system("cd $dir; cp -sf ../model_params.in .");
  } 

  # link multiphysics.in file
  if ( -e "multiphysics.in" ) {
    system("cd $dir; cp -sf ../multiphysics.in .");
  }
   
  print ("Finished creating data set for $nproc processors.\n\n");    
    
}


