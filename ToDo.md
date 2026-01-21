# ToDo Projekt

- [ ] Inflow und Outflow Boundary Conditions implementieren (Bene)
    - [ ] Inflow BC
    - [ ] Outflow BC
- [ ] Function als Geschwindigkeitsinput implementieren (für debugging purposes)

- [x] Precice integration testen (koppeln mit solverdummy) (julian)
  - [x] Dafür wrapper oder ein bisschen code architektur überlegen

- [ ] Herleitung der Kraftvektoren auf Festkörper (fabio)
    - Implementierung soll nur in eine Richtung kraft auswirken (orthogonal zur Wand, parallel vernachlässigen wir)
    - Was wenn Festkörper Rand berührt?

- [ ] Check Acceleration of precice to use in the advance method

- [x] Konfiguration überlegen von Tube
- [ ] Überlegung wie wir Zellen klassifizieren (fluid, solid, interface)

- [x] Calculix installieren
  - [x] Calculix Config Files checken und 2d tube erstellen (dafür bedarf es noch keiner preCICE koppellung) (julian)
  - [ ] Entrypoint für Calculix Simulation Schreiben (julian)

- [ ] Fenics installieren und testen

- [ ] Fluid-Solver und Structure Solver koppeln mit preCICE

- [ ] Plotting irgendwie (dicke von tube an stelle etc.), da mehr überlegen

- [ ] Präsentation vorbereiten
  - Schöne Bilder/Videos von Simulationen
  - Precice
  - Architektur Adapter
  - Herleitung Kraftvektoren
