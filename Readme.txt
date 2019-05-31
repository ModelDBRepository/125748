This entry to the ModelDB is associated to the paper entitled
"Diminished activity-dependent BDNF expression underlies cortical
neuron microcircuit hypoconnectivity resulting from exposure to mutant
huntingtin fragments", authored by Luca Gambazzi, Ozgun Gokce, Tamara
Seredenina, Elena Katsyuba, Heike Runne, Henry Markram, Michele
Giugliano and Ruth Luthi-Carter.

The model has been developed and implemented by Michele Giugliano
(mgiugliano@gmail.com)


Abstract:

Although previous studies of Huntington`s disease (HD) have addressed
many potential mechanisms of striatal neuron dysfunction and death, it
is also known based on clinical findings that cortical function is
dramatically disrupted in HD. With respect to disease etiology,
however, the specific molecular and neuronal circuit bases for the
cortical effects of mutant huntingtin (htt) have remained largely
unknown.

In the present work we studied the relation between the molecular
effects of mutant htt fragments in cortical cells and the
corresponding behavior of cortical neuron microcircuits using a novel
cellular model of HD. We observed that a transcript-selective
diminution in activity-dependent BDNF expression preceded the onset of
a synaptic connectivity deficit in ex vivo cortical networks, which
manifested as decreased spontaneous collective burst-firing behavior
measured by multi-electrode array substrates. Decreased BDNF
expression was determined to be a significant contributor to
network-level dysfunction, as shown by the ability of exogenous BDNF
to ameliorate cortical microcircuit burst firing.

The molecular determinants of the dysregulation of activity-dependent
BDNF expression by mutant htt appear to be distinct from previously
elucidated mechanisms, as they do not involve known
NRSF/REST-regulated promoter sequences, but instead result from
dysregulation of BDNF exon IV and VI transcription. These data
elucidate a novel HD-related deficit in BDNF gene regulation as a
plausible mechanism of cortical neuron hypoconnectivity and cortical
function deficits in HD. Moreover, the novel model paradigm
established here is well-suited to further mechanistic and drug
screening research applications.

A simple mathematical model is proposed to interpret the observations
and to explore the impact of specific synaptic dysfunctions on network
activity. Interestingly, the model predicts a decrease in synaptic
connectivity to be an early effect of mutant huntingtin in cortical
neurons, supporting the hypothesis of decreased, rather than
increased, synchronized cortical firing in HD.


Model information:

The model provided here is an ANSI-C source code, intended to be
compiled with "gcc" (native under Linux, by CygWin under Windows, by
Xcode under Mac OsX).

The source code is fully commented and rather simple to grasp, with
the aim of providing a useful reference to the interested user. The
software is intended (when compiled) to be invoked by some scripting
language (e.g. Matlab, or Octave). Libraries and includes are
provided. Matlab scripts for data analysis are also provided and
commented.


The simulation is launched by (Matlab) scripts (analysis*.m), which
provide a means to plot results as wells.


----------------------------------------------------------------------------

How to compile the code (from a 'shell' / 'terminal window' / et
similia):

(see also the source code at /giugliano/source/meanfield.c)


Let's assume the present working directory is 'giugliano (i.e. > pwd,
returns ..../giugliano)

compile by > " gcc -o meanfield source/meanfield.c -lm -O"

----------------------------------------------------------------------------

How to launch the simulation (from a 'shell' / 'terminal window' / et
similia):


Let's assume the present working directory is 'giugliano (i.e. > pwd,
returns ..../giugliano)


invoking the program with no inputÂ "meanfield" returns the usage..


USAGE: ./meanfield T N C I mext sext Use


where

T is the life time of the simulation

N is the number of (excitatory) neurons

C is the probability of random pairwise connection

I is the mean synaptic efficacy

mext, sext are statistical parameters related to the spontaneous
synaptic release

Use is the the release probability of short-term (facilitating and)
depressing synapses


Please refer to the published Supplemental Methods and to references
therein, or contact the authors.

----------------------------------------------------------------------------

How to 'play' with the parameters and extract meaningful information:


see the way 'schedule_simulation_example.m' is coded..


Michele Giugliano, PhD
