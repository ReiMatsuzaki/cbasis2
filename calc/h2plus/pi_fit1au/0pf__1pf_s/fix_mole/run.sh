CBASIS=../../../../..
ONEEIG=${CBASIS}/one_eig/bin/fast/one_eig
TWOPOT=${CBASIS}/two_pot/bin/fast/two_pot
${ONEEIG} one_eig.in.json | tee one_eig.out
${TWOPOT} two_pot.in.json | tee two_pot.out

