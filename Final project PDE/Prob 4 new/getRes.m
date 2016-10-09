function Res = getRes(interfW,interfS,interfE,interfWold,interfSold,interfEold)

        NormGreen = norm(interfW-interfWold);
        NormBlue  = norm([(interfS-interfSold);(interfE-interfEold)']);
        Res = max([norm(interfW-interfWold),norm([(interfS-interfSold);(interfE-interfEold)'])]);
end