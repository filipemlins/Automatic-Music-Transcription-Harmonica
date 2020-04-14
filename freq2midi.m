function midi = freq2midi(hz)

midi = 12*log2(hz/440)+69;