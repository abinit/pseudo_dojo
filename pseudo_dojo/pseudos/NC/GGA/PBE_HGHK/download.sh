root=http://cvs.berlios.de/cgi-bin/viewvc.cgi/cp2k/potentials/Goedecker/abinit/pbe

for p in `cat list.txt | cut -f 1`;
  do
  wget "${root}/${p}" -O ${p}.hgh.k
done

