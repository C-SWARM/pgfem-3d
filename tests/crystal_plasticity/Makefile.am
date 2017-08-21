noinst_SCRIPTS = 8.cpu.compression.sh \
                 8.cpu.tension.sh

TESTS = 8.cpu.compression.sh \
        8.cpu.tension.sh

8.cpu.compression.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "./common.sh 8 compression" >> $@
	chmod +x $@

8.cpu.tension.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "./common.sh 8 tension" >> $@
	chmod +x $@
