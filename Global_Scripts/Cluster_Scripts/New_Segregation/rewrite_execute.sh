#!/bin/bash
# This script rewrites the execute.sh file present in each architecture subfolder such that the seeds remain intact, but the commands after that are changed
# The commands that should be printed after the seeds are in execute_commands.txt:
commands_file=execute_commands.txt

for architecture in ./Arc*/
do
	# architecture=Arc1
	cd ${architecture}
	cp execute.sh old_execute.sh
	rm execute.sh
	# Preserving the seeds (first four lines) from the old execute script:
	cat > execute.sh <<- END
	$(head old_execute.sh -n 4)
	END
	cat "../${commands_file}" >> execute.sh
	rm old_execute.sh
	cd ..
done
