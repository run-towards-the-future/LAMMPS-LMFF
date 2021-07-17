#generate pc options
cpp -E $1 | grep -E 'PC[0-9]+_.*=\s*0x[0-9a-z]+' -o | awk -F '=' '{print "integer, parameter:: " $1 "=" $2}' | sed -e 's/0x\([0-9a-f]*\)/z"\1"/g' > $2
NPCR=$(grep -o -E "define\s+NPCR\s+[0-9]+" $1 | awk '{print $3}')
cat >> $2 <<EOF
type :: lwpf_evt_conf
integer(kind=8) :: pc_mask
EOF
for ((i=0;i<${NPCR};i++)); do
  echo "integer(kind=8):: evt_$i" >> $2
done
echo end type lwpf_evt_conf >> $2
