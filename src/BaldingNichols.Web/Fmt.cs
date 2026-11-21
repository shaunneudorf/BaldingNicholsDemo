using System.Globalization;

namespace BaldingNichols.Web;

/// <summary>Number formatting shared by the page and the formula components.</summary>
public static class Fmt
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Compact decimal for inputs to a formula, e.g. 0.1616 or 0.000304.</summary>
    public static string Plain(double v) => v.ToString("0.########", Inv);

    /// <summary>Six significant figures, switching to × 10ⁿ for very small values.</summary>
    public static string Prob(double v)
    {
        if (v == 0) return "0";
        if (Math.Abs(v) >= 1e-4) return v.ToString("0.######", Inv);
        return Sci(v);
    }

    /// <summary>Likelihood ratio: two decimals, or × 10ⁿ when very large.</summary>
    public static string Lr(double v) =>
        double.IsInfinity(v) || double.IsNaN(v) ? "undefined"
        : v < 1e6 ? v.ToString("N2", Inv)
        : Sci(v);

    private static string Sci(double v)
    {
        var s = v.ToString("0.####e+0", Inv);          // e.g. 3.0396e-5
        var parts = s.Split('e');
        var exp = int.Parse(parts[1], Inv);
        return $"{parts[0]} × 10{Superscript(exp)}";
    }

    private static string Superscript(int n)
    {
        const string map = "⁰¹²³⁴⁵⁶⁷⁸⁹";
        var chars = Math.Abs(n).ToString(Inv).Select(c => map[c - '0']);
        return (n < 0 ? "⁻" : "") + new string(chars.ToArray());
    }
}
