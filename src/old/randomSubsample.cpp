

virtual boost::shared_ptr<Population> randomSubsample(size_t size) const = 0;


boost::shared_ptr<Population> Population_Organisms::randomSubsample(size_t size) const
{
    if (size > organisms_.size())
        throw runtime_error("[Population::randomSubsample] Sample size exceeds population size.");

    set<size_t> indices;
    while (indices.size() < size) // may take a long time if size is close to organisms_.size()
        indices.insert(Random::uniform_integer(0,organisms_.size()-1));

    boost::shared_ptr<Population_Organisms> subsample(new Population_Organisms()); // TODO: fix

    for (set<size_t>::const_iterator it=indices.begin(); it!=indices.end(); ++it)
        subsample->organisms_.push_back(organisms_[*it]);

    return subsample;
}


