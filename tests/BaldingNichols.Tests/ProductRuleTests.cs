using BaldingNichols.Core;
using Xunit;

namespace BaldingNichols.Tests;

public class ProductRuleTests
{
    [Fact]
    public void Heterozygote_Is_2pq()
    {
        var r = BnCalculator.ProductRule(["12", "13"], new Dictionary<string, double> { ["12"] = 0.1, ["13"] = 0.2 });
        Assert.False(r.IsHomozygous);
        Assert.Equal(0.04, r.Probability, 12);
        Assert.Equal(25, r.LikelihoodRatio, 9);
    }

    [Fact]
    public void Homozygote_Is_PSquared()
    {
        var r = BnCalculator.ProductRule(["12", "12"], new Dictionary<string, double> { ["12"] = 0.1 });
        Assert.True(r.IsHomozygous);
        Assert.Equal(0.01, r.Probability, 12);
    }

    [Fact]
    public void NumericallyEquivalentAlleleNames_AreTheSameAllele()
    {
        var r = BnCalculator.Compute(0.01, ["09.3", "10"], ["9.3", "11"],
            new Dictionary<string, double> { ["9.3"] = 0.3, ["11"] = 0.2 });
        Assert.Equal(2, r.Orders[0].Draws[0].X);   // 9.3 seen once as victim, once as suspect
    }
}
