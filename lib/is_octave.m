function [isOctave]=is_octave

  isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

end
