#!/bin/bash
cp docs/usage.template usage.template
#looper --help > USAGE.temp 2>&1

for cmd in "--help" "init --help" "list --help" "listr --help" "pull --help" "build --help" "seek --help" "add --help" "remove --help" "getseq --help" "tag --help" "id --help" "subscribe --help" "unsubscribe --help"; do
	echo $cmd
	echo -e "## \`refgenie $cmd\`" > USAGE_header.temp
	refgenie $cmd --help > USAGE.temp 2>&1
	# sed -i 's/^/\t/' USAGE.temp
	sed -i.bak '1s;^;\`\`\`console\
;' USAGE.temp
#	sed -i '1s/^/\n\`\`\`console\n/' USAGE.temp
	echo -e "\`\`\`\n" >> USAGE.temp
	#sed -i -e "/\`looper $cmd\`/r USAGE.temp" -e '$G' usage.template  # for -in place inserts
	cat USAGE_header.temp USAGE.temp >> usage.template # to append to the end
done
rm USAGE.temp
rm USAGE_header.temp
rm USAGE.temp.bak
mv usage.template  docs/usage.md
cat docs/usage.md
#rm USAGE.temp
