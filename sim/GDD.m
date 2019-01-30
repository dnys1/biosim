function gdd = GDD(trackTa)
   % Calculate growing degree days(GDD)
   % https://en.wikipedia.org/wiki/Growing_degree-day
   Tbase = 10;
   T = trackTa - Tbase;
   T = T(T>0);
   gdd = max(trapz(T) / 24, 0);
end