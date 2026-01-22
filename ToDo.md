# ToDo Projekt

- [ ] Inflow und Outflow Boundary Conditions implementieren (Bene)
    - [x] Inflow BC
    - [x] Outflow BC
    - [ ] Testen ob das tut mit openFOAM

- [x] Function als Geschwindigkeitsinput implementieren (für debugging purposes)

- [x] Precice integration testen (koppeln mit solverdummy) (julian)
  - [x] Dafür wrapper oder ein bisschen code architektur überlegen

- [ ] Herleitung der Kraftvektoren auf Festkörper (fabio)
    - Implementierung soll nur in eine Richtung kraft auswirken (orthogonal zur Wand, parallel vernachlässigen wir)
    - Was wenn Festkörper Rand berührt?

- [ ] Check Acceleration of precice to use in the advance method

- [x] Konfiguration überlegen von Tube
- [ ] Überlegung wie wir Zellen klassifizieren (fluid, solid, interface)

- Calculix installieren
  - [x] julian
  - [ ] fabio
  - [ ] bene

- openFOAM installieren
  - [x] julian
  - [ ] fabio
  - [ ] bene

- [x] Calculix Config Files checken und 2d tube erstellen (dafür bedarf es noch keiner preCICE koppellung) (julian)
- [x] Entrypoint für Calculix Simulation Schreiben (julian)
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
