#!/bin/bash

function evo
{
#	dir=evodata2/$3.$1.`hexdump -n 4 /dev/urandom -e '"%08x"'`
	dir=evodata3/$3
	ssh -f $1 "cd src/evolution && mkdir -p $dir && nice -n19 ./src/crossover2 $dir evomany 5000 100 $2 .001 .2 .02 $3 </dev/null &>>$dir/log & disown"
	
}

# method xprob

evo asimov	 1 0
evo b2	 1 0.05
evo beckett	 1 0.1
evo dolly	 1 0.2
evo einstein	 1 0.3
evo fakiren	 1 0.4
evo glauser	 1 0.5
evo hank	 1 0.6
evo johnny	 1 0.7
evo kafka	 1 0.8
evo lev	 1 0.05
evo lindgren	 1 0.1
evo nabokov	 1 0.2
evo napoleon	 1 0.3
evo patsy	 1 0.4
evo piraten	 1 0.5
evo predator	 1 0.6
evo roth	 1 0.7
evo roy	 1 0.8
evo slas	 1 0
evo snorre	 1 0.05
evo st06	 1 0.1
evo tanya	 1 0.2
evo vilhelm	 1 0.3
