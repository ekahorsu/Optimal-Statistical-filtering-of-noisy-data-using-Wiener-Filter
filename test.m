filename = 'one_new.m4a';
[voice, Fs] = audioread(filename);
sound(voice, Fs)