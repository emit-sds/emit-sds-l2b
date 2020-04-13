# Warning, this is only to handle tetracorder's inability to read long lines
while read -r line; do sbatch -N 1 -n 1 -c 1 --mem=15000 --job-name=tetra_${line} -o ../salton_sea/l2b/logs/o_${line} -e ../salton_sea/l2b/logs/e_${line} -w node032 --wrap="sh run_tetracorder.sh ../salton_sea/l2a/${line}/output/${line}_rfl ../salton_sea/l2b/${line}_rfl_scaled ../${line}_rfl_scaled ../salton_sea/l2b/tetra_${line} ../salton_sea/l2b/${line}_sa"; done < ../salton_sea/lines.txt

# This is how it should be run
#while read -r line; do sbatch -N 1 -n 1 -c 20 --mem=50000 --job-name=tetra_${line} -o ../salton_sea/l2b/logs/o_${line} -e ../salton_sea/l2b/logs/e_${line} -w node032 --wrap="sh run_tetracorder.sh ${PWD}/../salton_sea/l2a/${line}/output/${line}_rfl ${PWD}/../salton_sea/l2b/${line}_rfl_scaled ${PWD}/../salton_sea/l2b/tetra_${line}"; done < ../salton_sea/lines.txt
