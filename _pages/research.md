---
permalink: /research/
title: "Research Topics"
author_profile: true
---

Thalamocortical circuit
======

Information from the sensory surround transmitts through the nuclei of the thalamus _en route_ to the cortex. 
The thalmic reticular nucleus sits between these structures and shapes thalamic activity through the inhibitory synapses sent to thalamus.

<iframe src="/files/images/TRN_circuit.png" width="100%" height="500" frameborder="no" border="0" marginwidth="0" marginheight="0"></iframe>
_credit: Julie Haas_

Electrical synapses
======

Electrical synapses are a major class of synapse formed by gap junctions between neurons.
These gap junctions allow for electrical current to flow between two coupled neurons, allowing voltage differences to influence the activity of the cells.


A simple model of an electrical synapse can thus be a static resistance applied to the voltage difference between two neurons:


$$
I_{cell1->cell2} = G_{elec}\cdot(V_{cell1}-V_{cell2})
$$

Computational models
======

My expertise with computational modelling is with utilizing the Hodgkin-Huxley (HH) based models.
Neuron spiking is simulated from the biophysics of the sodium and potassium ion channels and their gating properties, mathematically described by Hodgkin and Huxley. 

The simplest example is shown here:

<iframe src="/files/images/HHmodel.jpg" width="100%" height="500" frameborder="no" border="0" marginwidth="0" marginheight="0"></iframe>

A basic solving interface for this system could be written as:

```julia
using OrdinaryDiffEq

prob = ODEProblem(dsim!, u0, (startTime, endTime), p)
# dsim: HH function to be solved
# u0: initial conditions 
# p: parameters (channels, synapses etc.)

sol = solve(prob, BS3(), saveat=dt)
```

where parameters (p) are:

```julia
struct HHmodel
    gNa::Float64
    gK::Float64
    gL::Float64

    ENa::Float64
    EK::Float64
    EL::Float64

    C::Float64

    Iapp::Float64
end
```

and with the function `dsim!` describing the HH equations:

```julia
function dsim!(du, u, p, t)
	v, n, m, h = u

    I = (p.gK * (n^4.0) * (v-p.EK))
      + (p.gNa * (m^3.0) * h * (v-p.ENa))
      + (p.gL * (v-p.EL)) 
      + p.Iapp

    du[1]  = (-1.0/p.C) * I
    du[2]  = dn(v,n)
    du[3]  = dm(v,m)
    du[4]  = dh(v,h)
end
#dn,dm & dh dynamics ommited for clarity 
```

Solving the system with a positive square current illicits spiking such as:




Coming soon:
An interactive parameter explorer for neuron models
