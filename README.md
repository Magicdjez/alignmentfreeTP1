
Djeser Kordon
Wei Liu
# Résultats
Nous avons fait l'analyse des resultats des deux methodes dans le fichier analyser resultat. 
## Exec time / RAM

k = 21

### dictionnaires

```
/usr/bin/time -v python3 main.py 
index
GCF_000006945.2_ASM694v2_genomic.fna GCF_008244785.1_ASM824478v1_genomic.fna 0.9377423900044074 (0.9599864198488329, 0.9758858885062132)

	Command being timed: "python3 main.py"
	User time (seconds): 9.72
	System time (seconds): 0.29
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:10.02
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 934072
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 308576
	Voluntary context switches: 1
	Involuntary context switches: 84
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0


```

### Listes

```
/usr/bin/time -v python3 main.py 
lists
GCF_000006945.2_ASM694v2_genomic.fna GCF_008244785.1_ASM824478v1_genomic.fna 0.9377423900044074 (0.9599864198488329, 0.9758858885062132)

	Command being timed: "python3 main.py"
	User time (seconds): 9.04
	System time (seconds): 0.16
	Percent of CPU this job got: 100%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:09.20
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 567176
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 163235
	Voluntary context switches: 1
	Involuntary context switches: 70
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

```

