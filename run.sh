for step in 2 4 8 16 32 64 128 256 512 1024
do
#  for order in 4 6 8 10 12
  for order in 14
  do
    for npec in 1 2 3
#	rm plot/${order}th-pec${npec}.dat
	do
      ./eight-${order}th ${step} ${npec} | grep grep >> plot/${order}th-pec${npec}.dat
	done
  done
#./eight-12th $step 1 | grep grep >> 12th-pec1.dat
#  ./eight-12th $step 2 | grep grep >> 12th-pec2.dat
#  ./eight-12th $step 3 | grep grep >> 12th-pec3.dat
done
