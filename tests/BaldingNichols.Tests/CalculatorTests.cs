using BaldingNichols.Core;
using Xunit;

namespace BaldingNichols.Tests;

public class CalculatorTests
{
    private static Dictionary<string, double> Freqs(params (string, double)[] f) =>
        f.ToDictionary(x => x.Item1, x => x.Item2);

    [Fact]
    public void SingleSource_Heterozygote_MatchesNrcII_4_2()
    {
        double t = 0.03, p = 0.1, q = 0.2;
        var r = BnCalculator.Compute(t, [], ["12", "13"], Freqs(("12", p), ("13", q)));
        double expected = 2 * (t + (1 - t) * p) * (t + (1 - t) * q) / ((1 + t) * (1 + 2 * t));
        Assert.Equal(expected, r.GenotypeProbability, 12);
    }

    [Fact]
    public void SingleSource_Homozygote_MatchesNrcII_4_2()
    {
        double t = 0.03, p = 0.1;
        var r = BnCalculator.Compute(t, [], ["12", "12"], Freqs(("12", p)));
        double expected = (2 * t + (1 - t) * p) * (3 * t + (1 - t) * p) / ((1 + t) * (1 + 2 * t));
        Assert.Equal(expected, r.GenotypeProbability, 12);
        Assert.Single(r.Orders);
    }

    [Fact]
    public void ThetaZero_GivesProductRule_EvenWithVictim()
    {
        var f = Freqs(("12", 0.15), ("13", 0.25));
        var r = BnCalculator.Compute(0, ["11", "12"], ["12", "13"], f);
        Assert.Equal(2 * 0.15 * 0.25, r.GenotypeProbability, 12);
        Assert.Equal(BnCalculator.ProductRuleProbability(["12", "13"], f), r.GenotypeProbability, 12);
    }

    [Fact]
    public void TwoPersonExample_MatchesSpreadsheet()
    {
        double p = 0.161398176291793;
        var r = BnCalculator.Compute(0.01, ["11", "12"], ["12", "13"], Freqs(("12", p), ("13", p)));
        Assert.Equal(0.05699125216039515, r.GenotypeProbability, 12);
        Assert.Equal(17.546552533810242, r.LikelihoodRatio, 9);
        Assert.Equal(2, r.Orders.Count);
        Assert.Equal(r.Orders[0].Probability, r.Orders[1].Probability, 12);
    }

    [Fact]
    public void Homozygote_WithVictim_UsesOneOrderAndSingleFrequency()
    {
        double t = 0.01, p = 0.2;
        var r = BnCalculator.Compute(t, ["11", "12"], ["12", "12"], Freqs(("12", p)));
        double a = (3 * t + (1 - t) * p) / (1 + 3 * t);   // x=3, n=4
        double b = (4 * t + (1 - t) * p) / (1 + 4 * t);   // x=4, n=5
        Assert.Equal(a * b, r.GenotypeProbability, 12);
        Assert.Single(r.Orders);
    }

    [Fact]
    public void BlankVictim_IsTreatedAsNoVictim()
    {
        var f = Freqs(("12", 0.1), ("13", 0.2));
        var a = BnCalculator.Compute(0.02, ["", ""], ["12", "13"], f);
        var b = BnCalculator.Compute(0.02, [], ["12", "13"], f);
        Assert.Equal(b.GenotypeProbability, a.GenotypeProbability, 15);
    }

    [Theory]
    [InlineData(-0.1)]
    [InlineData(1.0)]
    public void InvalidTheta_Throws(double theta)
    {
        Assert.Throws<ArgumentException>(() =>
            BnCalculator.Compute(theta, [], ["12", "13"], Freqs(("12", 0.1), ("13", 0.2))));
    }

    [Fact]
    public void MissingFrequency_Throws()
    {
        Assert.Throws<ArgumentException>(() =>
            BnCalculator.Compute(0.01, [], ["12", "13"], Freqs(("12", 0.1))));
    }
}
