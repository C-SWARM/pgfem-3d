This is a simple unit cell containing a stiff periodic particle in a
softer binder.

### Modifying the geometry ###
If you chose to modify the geometry in `pack.rop2t3d`, you will need
to manually run `RoP2T3d` for them to take effect. Furthermore, you
will need to ensure that the corner node ids and region ids have not
changed, or update `pack.out.header` accordingly.

### Generating input deck ###
To generate the input deck, do:
```
./makeset.pl -np <number of processors>
```
To clean (delete) all files created by the `makeset.pl` script, do:
```
./makest.pl clean
```

### Running the example ###
Modify the `runit.sh` script to point to the appropriate version of
the code (default: master). Then build the input deck for the desired
number of cores (see previous section). Finally, do:
```
./runit.sh <number of processors>
```
