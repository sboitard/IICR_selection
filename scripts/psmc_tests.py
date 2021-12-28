import subprocess

reps = 10
L = int(1e8)
mu = 2.5e-8
N = 1000

N_classes = {
    (1, 10) : 'A',
    # (0.1, 1) : 'B'
}
N_proportions = {
    (0.50, 0.50) : 'A',
    # (0.90, 0.10) : 'B',
    # (0.99, 0.01) : 'C'
}
chunks = {
    # int(1e6) : 'A',
    # int(1e5) : 'B',
    # int(1e4) : 'C'
    int(5e3) : 'D',
    int(2e3) : 'E',
    int(1e3) : 'F'
}

psmc_commands_file = open(f"./psmcs_to_run.txt", "w")

for (n_class, n_class_tag) in N_classes.items():
    for (n_proportion, proportion_tag) in N_proportions.items():
        for (chunk, chunk_tag) in chunks.items():
            tag = f'_s{n_class_tag}_f{proportion_tag}_c{chunk_tag}'

            for i in range(len(n_class)):
                Ne = int(n_class[i] * N)
                Le = int(n_proportion[i] * L)
                theta = 4 * mu * Le * Ne
                rho = theta / 5

                ms_command = f"./ms 2 {reps} -t {theta} -r {rho} {Le} -p 8"
                ms_output = open(f'./tests_class_{i}{tag}.ms', 'w')
                print("running", ms_command)
                subprocess.run(ms_command.split(), check = True, stdout = ms_output)
                
                print("running ms2psmcfa")
                psmcfa_output = open(f'./tests_class_{i}{tag}.psmcfa', 'w')
                subprocess.run(f"python3 ms2psmcfa.py tests_class_{i}{tag}.ms".split(), check = True, stdout = psmcfa_output)
                print("done!")

            bins = 100
            fa_files = [
                f"./tests_class_0{tag}.psmcfa",
                f"./tests_class_1{tag}.psmcfa",
            ]

            fa_data = []

            for fa_file in fa_files:
                current = open(fa_file, "r").readlines()
                current = ''.join(current)
                current = current.replace('\n', '')
                current = current.split('>')
                fa_data.append(current)

            final_data = []

            # compressed chunks
            cc = [ 
                int((n_proportion[0] * chunk / n_proportion[1]) / bins),
                int(chunk / bins)
            ]
            parts = int(L * n_proportion[1] / chunk)

            for i in range(reps):
                class_a = fa_data[0][i + 1][1:]
                class_b = fa_data[1][i + 1][1:]
                final_data.append(f">{i + 1}\n")

                chromosome = ''
                for j in range(parts):
                    chromosome += class_a[cc[0] * j : cc[0] * (j + 1)]
                    chromosome += class_b[cc[1] * j : cc[1] * (j + 1)]
                    
                final_data.append(chromosome + "\n")

            final_file = open(f"./results{tag}.psmcfa", "w")
            final_file.writelines(final_data)
            psmc_command = f'./psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o results{tag}.psmc results{tag}.psmcfa'
            psmc_commands_file.writelines([psmc_command + '\n'])
            print("PSMC command exported\n\n")
            # print("PSMC command:")
            # print(f'./psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o results{tag}.psmc results{tag}.psmcfa')
            # subprocess.run(f'./psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o results{tag}.psmc results{tag}.psmcfa'.split(), check=True)