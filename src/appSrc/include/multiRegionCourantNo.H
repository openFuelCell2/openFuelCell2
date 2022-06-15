    scalar CoNum = -GREAT;

    forAll(regions, i)
    {
        if(regions[i].type() == "fluid")
        {
            const phaseSystem& phase = regions[i].lookupObject<phaseSystem>(phaseSystem::propertiesName);

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
