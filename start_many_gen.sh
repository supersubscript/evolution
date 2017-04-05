#!/bin/bash



function evo
{
	dir=evodata2/
	ssh -f $1 "cd src/evolution && { \
	nice -n19 ./src/crossover2 $dir generations -30 2 .001 .2 .02 $2; \
	}  "
#	}  </dev/null &>>$dir/log & disown"
}


exit
evo lindgren	0
evo merle		0.05
evo johnny		0.1
evo asimov		0.2
evo beckett		0.3
evo predator	0.4
evo roy			0.5
evo tanya		0.6

