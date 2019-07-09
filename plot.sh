 #!/bin/bash
conftypes=${1:-s064t16_b0687361 s080t20_b0703500 s096t24_b0719200}
qcdtype=${2:-quenched}


for conftype in $conftypes ; do
    beta=${conftype#*_b}; beta=`bc <<< "scale=5;$beta/100000"`
    ns=${conftype#s}; ns=${ns%%t*}
    nt=${conftype#*t}; nt=${nt%%_b*}
    python3 plot.py $qcdtype $conftype $beta $ns $nt
done
