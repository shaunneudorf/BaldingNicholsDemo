using BaldingNichols.Core;
using Xunit;

namespace BaldingNichols.Tests;

public class AlleleDatabaseTests
{
    // L1: alleles 8, 9, 9.3 observed (k = 3, N = 100).  L2: allele 8 not observed (k = 2, N = 50).
    private const string Csv = "Allele,L1,L2\r\n8,0.5,0\r\n9,0.25,0.4\r\n9.3,0.25,0.6\r\nN,100,50\r\n";

    private static AlleleDatabase Db() => AlleleDatabase.ParseCsv(Csv);

    [Fact]
    public void Parse_ReadsLociNAndK()
    {
        var db = Db();
        Assert.Equal(new[] { "L1", "L2" }, db.Loci.Select(l => l.Name).ToArray());
        Assert.Equal(100, db.Find("L1")!.N);
        Assert.Equal(3, db.Find("L1")!.K);
        Assert.Equal(50, db.Find("L2")!.N);
        Assert.Equal(2, db.Find("l2")!.K);
    }

    // ---- Balding-Nichols frequency: (x + 1/k) / (N + 1) ----

    [Fact]
    public void BaldingNichols_ObservedAllele_UsesCountPlusOneOverK()
    {
        // L1: N = 100, k = 3, allele 9.3 has frequency 0.25 so x = 25.
        var r = Db().Find("L1")!.LookupBaldingNichols("9.3");
        Assert.True(r.Observed);
        Assert.Equal(25, r.Count);
        Assert.Equal((25 + 1.0 / 3) / 101, r.Value, 12);
    }

    [Fact]
    public void BaldingNichols_UnseenAllele_IsTheXEqualsZeroCase()
    {
        var l1 = Db().Find("L1")!;
        var r = l1.LookupBaldingNichols("40");
        Assert.False(r.Observed);
        Assert.Equal((1.0 / 3) / 101, r.Value, 12);
        Assert.Equal(l1.UnseenBaldingNicholsFrequency, r.Value, 15);
    }

    [Fact]
    public void BaldingNichols_ZeroFrequencyRowIsUnseen()
    {
        var r = Db().Find("L2")!.LookupBaldingNichols("8");
        Assert.False(r.Observed);
        Assert.Equal((1.0 / 2) / 51, r.Value, 12);
    }

    // ---- Product rule: database frequency, or 5/2N when unseen ----

    [Fact]
    public void ProductRule_ObservedAllele_IsTheRawDatabaseFrequency()
    {
        var r = Db().Find("L1")!.LookupProductRule("9.3");
        Assert.True(r.Observed);
        Assert.Equal(0.25, r.Value, 12);
    }

    [Fact]
    public void ProductRule_UnseenAllele_UsesFiveOverTwoN()
    {
        // N is the allele count (2N individuals' worth), so 5/2N = 5/100 for L1 and 5/50 for L2.
        Assert.Equal(5.0 / 100, Db().Find("L1")!.LookupProductRule("40").Value, 12);
        Assert.Equal(5.0 / 50, Db().Find("L2")!.LookupProductRule("8").Value, 12);
        Assert.False(Db().Find("L2")!.LookupProductRule("8").Observed);
    }

    [Theory]
    [InlineData("09.3")]
    [InlineData("9.30")]
    [InlineData(" 9.3 ")]
    public void Lookup_MatchesNumerically(string typed)
    {
        Assert.True(Db().Find("L1")!.LookupBaldingNichols(typed).Observed);
        Assert.True(Db().Find("L1")!.LookupProductRule(typed).Observed);
    }

    [Fact]
    public void MissingNRow_Throws()
    {
        Assert.Throws<FormatException>(() => AlleleDatabase.ParseCsv("Allele,L1\r\n8,1\r\n"));
    }

    [Fact]
    public void RealDatabase_LoadsAndEachLocusSumsToOne()
    {
        var db = LoadReal();

        Assert.Equal(15, db.Loci.Count);

        var d8 = db.Find("D8S1179")!;
        Assert.Equal(328, d8.N);
        Assert.Equal(10, d8.K);
        Assert.Equal(326, db.Find("vWA")!.N);

        foreach (var l in db.Loci)
            Assert.InRange(l.Observed.Sum(kv => kv.Value), 0.999, 1.001);
    }

    [Fact]
    public void RealDatabase_FrequenciesConvertToWholeCountsThatSumToN()
    {
        foreach (var l in LoadReal().Loci)
        {
            foreach (var kv in l.Observed)
                Assert.InRange(Math.Abs(kv.Value * l.N - l.Count(kv.Key)), 0, 0.05);
            Assert.Equal(l.N, l.Observed.Sum(kv => l.Count(kv.Key)));
        }
    }

    [Fact]
    public void RealDatabase_BaldingNicholsFrequenciesOfObservedAllelesSumToOne()
    {
        // sum of (x + 1/k)/(N+1) over the k observed alleles = (N + 1)/(N + 1)
        foreach (var l in LoadReal().Loci)
            Assert.Equal(1.0, l.Observed.Sum(kv => l.LookupBaldingNichols(kv.Key).Value), 9);
    }

    [Fact]
    public void UnseenAllele_FlowsThroughCalculator()
    {
        var l1 = Db().Find("L1")!;
        var freqs = new Dictionary<string, double>
        {
            ["9"] = l1.LookupBaldingNichols("9").Value,
            ["40"] = l1.LookupBaldingNichols("40").Value,
        };
        var r = BnCalculator.Compute(0.01, [], ["9", "40"], freqs);
        Assert.Equal((1.0 / 3) / 101, r.Orders[0].Draws[1].Frequency, 12);
    }

    private static AlleleDatabase LoadReal() =>
        AlleleDatabase.ParseCsv(File.ReadAllText(Path.Combine(AppContext.BaseDirectory, "Caucasian_database.csv")));
}
