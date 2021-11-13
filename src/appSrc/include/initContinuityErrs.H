HashTable<scalar> cumulativeContErr(fluidRegion.size(), 0.0);

forAllConstIter(wordList, fluid.toc(), iter)
{
    cumulativeContErr.insert
    (
        *iter,
        0.0
    );
}
