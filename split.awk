#!/usr/bin/gawk -f

BEGIN{RS="&\n"; FS="\n"}
{
  # print "FOUND FIELD: "$1
  if($1 ~ /CONTROL/) {
    for(i=2; i<=NF; ++i)
    {
      print $i >> "CTRL"
    }
  }
  else if($1 ~ /SYSTEM/) {
    for(i=2; i<=NF; ++i) {
      print $i >> "SYSTEM"
    }
  } else if ($1 ~ /ELECTRONS/) {
    for(i=2; i<=NF; ++i) {
      print $i >> "ELECTRONS"
    }
  } else if ($i ~ /SPECIES/){
    for(i=2; i<=NF; ++i) {
      print $i >> "SPECIES"
    }
  } else if ($1 ~ /ATOMIC_POSITIONS/) {
    for(i=2; i<=NF; ++i) {
      print $i >> "POS"
    }
  } else if ($1 ~ /ATOMIC_SPECIES/) {
    for(i=2; i<=NF; ++i) {
      print $i >> "SPECIES"
    }
  }else if ($1 ~ /K_POINTS/) {
    for (i=2;i<=NF;++i) {
      print $i >> "KPOINTS"
    }
  } else if ($1 ~ /CELL_PARAMETERS/) {
    for (i=2;i<=NF;++i) {
      print $i >> "CELL"
    }
  }
}
END{}
