# Seawater reservoir, with proton pumping and Ca influx

1D model

## Equations

```math
\frac{d[Ca]_{CF}}{dt} = \frac{0.5 J}{M_{CF}} + \frac{F_{SW}}{M_{CF}} ([Ca]_{SW} - [Ca]_{CF}) - \frac{G}{M_{CF}} \\
\frac{dTA_{CF}}{dt} = \frac{J}{M_{CF}} + \frac{F_{SW}}{M_{CF}} (TA_{SW} - TA_{CF}) - \frac{2G}{M_{CF}} + P_{HCO_3} \\
\frac{dDIC_{CF}}{dt} = \frac{F_{SW}}{M_{CF}} (DIC_{SW} - DIC_{CF}) - \frac{G}{M_{CF}} + D_{CO2} + P_{HCO_3} \\
G = k (\Omega - 1)^{n} \\
```

Where $G$ is calcification rate (mol cm-2 s-1), $J$ is proton pumping rate (mol cm-2 s-1), and $F_{SW}$ is the flux of seawater (kg cm-2 s-1).

Precipitation rate can be calculated as $k(\Omega - 1)^n$ where, for calculation speed, we calculate $[CO_3^{2-}]$ as `min((TA - DIC), DIC).`

```math
\frac{d[TE]_{CF}}{dt} = \frac{F_{SW}}{M_{CF}} ([TE]_{SW} - [TE]_{CF}) - \frac{K^{G}_{TE} ([TE]_{CF} / [Ca]_{CF}) G}{M_{CF}} + \frac{K^{J}_{TE} ([TE]_{SW} / [Ca]_{SW}) J}{M_{CF}}
```

Igoring DIC for simplicity, this reduces to:

```math
\frac{d[Ca]_{CF}}{dt} = \frac{0.5 J}{M_{CF}} + \frac{F_{SW}}{M_{CF}} ([Ca]_{SW} - [Ca]_{CF}) - \frac{G}{M_{CF}} \\
\frac{dTA_{CF}}{dt} = \frac{J}{M_{CF}} + \frac{F_{SW}}{M_{CF}} (TA_{SW} - TA_{CF}) - \frac{2G}{M_{CF}}\\
\frac{d[TE]_{CF}}{dt} = \frac{F_{SW}}{M_{CF}} ([TE]_{SW} - [TE]_{CF}) - \frac{K^{G}_{TE} ([TE]_{CF} / [Ca]_{CF}) G}{M_{CF}} + \frac{K^{J}_{TE} ([TE]_{SW} / [Ca]_{SW}) J}{M_{CF}}
```

Other sensitivies to include:

```math
K^{G}_{TE} = K^{G}_{ref} e^{k^G (T - T_{ref})}
```

where $K^{G}_{ref}$ is the inorganic partition coefficient at a reference temperature $T_{ref}$, and $k^G$ is the temperature sensitivity of the partition coefficient (0.03 for inorganic calcite).

## Alkalinity Counter-Ion Considerations

The above assumes that all protons are pumped via a Ca-ATPase pump with stoichiometry $H^+:0.5Ca^{2+}$. This may not be the case, and we can introduce a parameter, $f_{Ca}$, to account for this, where $f_{Ca}$ is the fraction of protons pumped via a Ca-ATPase pump. The remaining protons are pumped via a different pump that uses another ion (e.g. Na) to balance charge. This modifies the equation to:

```math
\frac{d[Ca]_{CF}}{dt} = \frac{0.5 J \color{red}{f_{Ca}}}{M_{CF}} + \frac{F_{SW}}{M_{CF}} ([Ca]_{SW} - [Ca]_{CF}) - \frac{G}{M_{CF}} \\
\frac{dTA_{CF}}{dt} = \frac{J}{M_{CF}} + \frac{F_{SW}}{M_{CF}} (TA_{SW} - TA_{CF}) - \frac{2G}{M_{CF}}\\
\frac{d[TE]_{CF}}{dt} = \frac{F_{SW}}{M_{CF}} ([TE]_{SW} - [TE]_{CF}) - \frac{K^{G}_{TE} ([TE]_{CF} / [Ca]_{CF}) G}{M_{CF}} + \frac{K^{J}_{TE} ([TE]_{SW} / [Ca]_{SW}) J}{M_{CF}}
```

This should allow the model to simulate calcification without a massive over-transport of Ca!