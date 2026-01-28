# ToDo Projekt

- [ ] Inflow und Outflow Boundary Conditions implementieren (Bene)
    - [x] Inflow BC
    - [x] Outflow BC
    - [x] Testen ob das tut mit openFOAM (julian)

- [x] Function als Geschwindigkeitsinput implementieren (für debugging purposes)

- [x] Precice integration testen (koppeln mit solverdummy) (julian)
  - [x] Dafür wrapper oder ein bisschen code architektur überlegen

- [ ] Herleitung der Kraftvektoren auf Festkörper (fabio)
    - Implementierung soll nur in eine Richtung kraft auswirken (orthogonal zur Wand, parallel vernachlässigen wir)
    - Was, wenn Festkörper Rand berührt?

- [ ] Check Acceleration of precice to use in the advance method

- [x] Konfiguration überlegen von Tube (julian)

- [ ] Überlegung wie wir Zellen klassifizieren (fluid, solid, interface) (bene)

- Calculix installieren
  - [x] julian
  - [ ] fabio
  - [x] bene

- openFOAM installieren
  - [x] julian
  - [ ] fabio
  - [ ] bene

- [x] Calculix Config Files checken und 2d tube erstellen (dafür bedarf es noch keiner preCICE koppellung) (julian)
- [x] Entrypoint für Calculix Simulation Schreiben (julian)
- [x] OpenFOAM entrypoint schreiben (would be nice to read our fluid config and simulate same scenario with openFOAM) (julian)
    - Wurde uns auch im Tutorium ans Herz gelegt um dann in der Präsentation unseren solver gegen openFOAM plotten zu können ähnlich zu den preCICE tutorials 
- [ ] Calculix Coupling testen mit OpenFOAM (julian)

- [ ] Fluid-Solver und Calculix koppeln mit preCICE
  - Hängt von Kraftvektoren und Boundary conditions ab
  - [ ] Simulation als API callbar machen

- [ ] Plotting irgendwie (dicke von tube an stelle etc.), da mehr überlegen

- [ ] Präsentation vorbereiten
  - Schöne Bilder/Videos von Simulationen
  - Precice
  - Architektur Adapter
  - Herleitung Kraftvektoren
