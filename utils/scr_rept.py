def prepare(max_repeat):
    global As, Cs, Gs, Ts
    As = '0' * (max_repeat+1)
    Cs = '1' * (max_repeat+1)
    Gs = '2' * (max_repeat+1)
    Ts = '3' * (max_repeat+1)

def screen_repeat(drop, max_repeat, gc_dev):
    dna = drop.toDNA()
    if As in dna or Cs in dna or Gs in dna or Ts in dna: 
        return 0
    gc = dna.count("1") + dna.count("2")  
    gc = gc/(len(dna)+0.0)
    if (gc < 0.5 - gc_dev) or (gc > 0.5 + gc_dev):
        return 0
    return 1