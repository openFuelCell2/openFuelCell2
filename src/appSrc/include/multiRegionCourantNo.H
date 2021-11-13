    scalar CoNum = -GREAT;

    forAll(regions, i)
    {
        if(regions[i].type() == "fluid")
        {
            const phaseSystem& phase = regions[i].lookupObject<phaseSystem>("phaseProperties");

            CoNum = max
            (
                regionCourantNo
                (
                    regions[i],
                    runTime,
                    phase.phirMag()
                ),
                CoNum
            );
        }
    }
