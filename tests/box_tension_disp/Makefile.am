check_SCRIPTS=1.cpu.sh \
              2.cpu.sh \
              4.cpu.sh \
              8.cpu.sh \
             16.cpu.sh 

TESTS=1.cpu.sh \
      2.cpu.sh \
      4.cpu.sh \
      8.cpu.sh \
     16.cpu.sh 

1.cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "./run.sh 1" >> $@
	chmod +x $@

2.cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "./run.sh 2" >> $@
	chmod +x $@

4.cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "./run.sh 4" >> $@
	chmod +x $@

8.cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "./run.sh 8" >> $@
	chmod +x $@

16.cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "./run.sh 16" >> $@
	chmod +x $@
