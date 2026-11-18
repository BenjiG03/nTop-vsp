using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;
using System.Text;
using Accessibility;

// Drives a WinForms window through MSAA and window messages, so nothing depends on which
// window has the keyboard. Only windows of the given process are ever touched.
public static class Msaa
{
    [DllImport("oleacc.dll")] static extern int AccessibleObjectFromWindow(IntPtr hwnd, uint id, ref Guid iid, [MarshalAs(UnmanagedType.IUnknown)] out object obj);
    [DllImport("oleacc.dll")] static extern int AccessibleChildren(IAccessible c, int start, int count, [Out, MarshalAs(UnmanagedType.LPArray, SizeParamIndex = 2)] object[] kids, out int got);
    [DllImport("user32.dll")] static extern bool EnumWindows(EnumProc cb, IntPtr l);
    [DllImport("user32.dll")] static extern bool EnumChildWindows(IntPtr parent, EnumProc cb, IntPtr l);
    [DllImport("user32.dll")] static extern uint GetWindowThreadProcessId(IntPtr h, out uint pid);
    [DllImport("user32.dll", CharSet = CharSet.Unicode)] static extern int GetWindowText(IntPtr h, StringBuilder s, int n);
    [DllImport("user32.dll", CharSet = CharSet.Unicode)] static extern int GetClassName(IntPtr h, StringBuilder s, int n);
    [DllImport("user32.dll")] static extern bool IsWindowVisible(IntPtr h);
    [DllImport("user32.dll")] static extern bool IsWindowEnabled(IntPtr h);
    [DllImport("user32.dll")] static extern int GetDlgCtrlID(IntPtr h);
    [DllImport("user32.dll")] static extern bool GetWindowRect(IntPtr h, out RECT r);
    [DllImport("user32.dll")] static extern bool PrintWindow(IntPtr h, IntPtr dc, uint flags);
    [DllImport("user32.dll", CharSet = CharSet.Unicode, EntryPoint = "SendMessageW")] static extern IntPtr SendText(IntPtr h, uint m, IntPtr w, string l);
    [DllImport("user32.dll", CharSet = CharSet.Unicode, EntryPoint = "SendMessageW")] static extern IntPtr SendBuf(IntPtr h, uint m, IntPtr w, StringBuilder l);
    [DllImport("user32.dll")] static extern IntPtr SendMessageTimeout(IntPtr h, uint m, IntPtr w, IntPtr l, uint flags, uint ms, out IntPtr result);
    [DllImport("user32.dll")] static extern bool PostMessage(IntPtr h, uint m, IntPtr w, IntPtr l);
    public delegate bool EnumProc(IntPtr h, IntPtr l);
    [StructLayout(LayoutKind.Sequential)] public struct RECT { public int L, T, R, B; }

    const uint OBJID_CLIENT = 0xFFFFFFFC, OBJID_WINDOW = 0;
    const uint WM_SETTEXT = 0x000C, WM_GETTEXT = 0x000D, WM_GETTEXTLENGTH = 0x000E, WM_CLOSE = 0x0010;
    const uint BM_CLICK = 0x00F5, CB_SETCURSEL = 0x014E, WM_LBUTTONDOWN = 0x0201, WM_LBUTTONUP = 0x0202;

    public static List<IntPtr> TopWindows(int pid)
    {
        var found = new List<IntPtr>();
        EnumWindows((h, l) =>
        {
            uint p; GetWindowThreadProcessId(h, out p);
            if (p == (uint)pid && IsWindowVisible(h)) found.Add(h);
            return true;
        }, IntPtr.Zero);
        return found;
    }

    public static string Text(IntPtr h)
    {
        IntPtr len;
        if (SendMessageTimeout(h, WM_GETTEXTLENGTH, IntPtr.Zero, IntPtr.Zero, 2, 2000, out len) == IntPtr.Zero) return "<busy>";
        var s = new StringBuilder(len.ToInt32() + 2);
        SendBuf(h, WM_GETTEXT, (IntPtr)s.Capacity, s);
        return s.ToString();
    }

    public static string Class(IntPtr h) { var s = new StringBuilder(256); GetClassName(h, s, 256); return s.ToString(); }

    public static string Describe(IntPtr top)
    {
        var o = new StringBuilder();
        RECT r; GetWindowRect(top, out r);
        o.AppendFormat("TOP {0} '{1}' class={2} [{3},{4},{5},{6}]\n", top, Text(top), Class(top), r.L, r.T, r.R - r.L, r.B - r.T);
        EnumChildWindows(top, (h, l) =>
        {
            RECT c; GetWindowRect(h, out c);
            string cls = Class(h);
            int i = cls.IndexOf(".app.");
            if (i > 0) cls = cls.Substring(0, i).Replace("WindowsForms10.", "");
            string t = Text(h);
            if (t.Length > 60) t = t.Substring(0, 60) + "...";
            o.AppendFormat("  {0} id={1} {2} vis={3} en={4} [{5},{6},{7},{8}] '{9}'\n", h, GetDlgCtrlID(h), cls, IsWindowVisible(h) ? 1 : 0, IsWindowEnabled(h) ? 1 : 0, c.L - r.L, c.T - r.T, c.R - c.L, c.B - c.T, t.Replace("\r", " ").Replace("\n", " "));
            return true;
        }, IntPtr.Zero);
        return o.ToString();
    }

    public static int CtrlId(IntPtr h) { return GetDlgCtrlID(h); }

    public static List<IntPtr> Children(IntPtr top)
    {
        var l2 = new List<IntPtr>();
        EnumChildWindows(top, (h, l) => { l2.Add(h); return true; }, IntPtr.Zero);
        return l2;
    }

    static IAccessible Acc(IntPtr h, uint id)
    {
        var iid = new Guid("618736E0-3C3D-11CF-810C-00AA00389B71");
        object o;
        int hr = AccessibleObjectFromWindow(h, id, ref iid, out o);
        if (hr != 0) throw new Exception("AccessibleObjectFromWindow " + hr);
        return (IAccessible)o;
    }

    static object[] Kids(IAccessible a)
    {
        int n = a.accChildCount;
        var kids = new object[n];
        int got = 0;
        if (n > 0) AccessibleChildren(a, 0, n, kids, out got);
        var r = new object[got];
        Array.Copy(kids, r, got);
        return r;
    }

    static void Walk(IAccessible a, object child, int depth, int maxDepth, StringBuilder o)
    {
        string name = "", value = ""; int role = 0, state = 0;
        try { name = a.get_accName(child); } catch { }
        try { value = a.get_accValue(child); } catch { }
        try { role = Convert.ToInt32(a.get_accRole(child)); } catch { }
        try { state = Convert.ToInt32(a.get_accState(child)); } catch { }
        o.AppendFormat("{0}role={1} state={2:x} name='{3}' value='{4}'\n", new string(' ', 2 * depth), role, state, name, value);
        if (depth >= maxDepth || !(child is int) || (int)child != 0) return;
        foreach (var k in Kids(a))
        {
            if (k is IAccessible) Walk((IAccessible)k, 0, depth + 1, maxDepth, o);
            else Walk(a, k, depth + 1, maxDepth, o);
        }
    }

    public static string DumpAcc(IntPtr h, int maxDepth, bool client)
    {
        var o = new StringBuilder();
        Walk(Acc(h, client ? OBJID_CLIENT : OBJID_WINDOW), 0, 0, maxDepth, o);
        return o.ToString();
    }

    // Finds the first accessible with this exact name under h and performs its default action.
    static bool Press(IAccessible a, object child, string want, int depth)
    {
        string name = null;
        try { name = a.get_accName(child); } catch { }
        if (name != null && name.Replace("&", "") == want) { a.accDoDefaultAction(child); return true; }
        if (depth > 10 || !(child is int) || (int)child != 0) return false;
        foreach (var k in Kids(a))
        {
            bool done = k is IAccessible ? Press((IAccessible)k, 0, want, depth + 1) : Press(a, k, want, depth + 1);
            if (done) return true;
        }
        return false;
    }

    public static bool PressName(IntPtr h, string want) { return Press(Acc(h, OBJID_WINDOW), 0, want, 0); }

    public static void SetText(IntPtr h, string s) { SendText(h, WM_SETTEXT, IntPtr.Zero, s); }
    public static void Click(IntPtr h) { PostMessage(h, BM_CLICK, IntPtr.Zero, IntPtr.Zero); }
    public static void Close(IntPtr h) { PostMessage(h, WM_CLOSE, IntPtr.Zero, IntPtr.Zero); }

    public static void Picture(IntPtr h, string path)
    {
        RECT r; GetWindowRect(h, out r);
        using (var bmp = new Bitmap(r.R - r.L, r.B - r.T))
        using (var g = Graphics.FromImage(bmp))
        {
            IntPtr dc = g.GetHdc();
            PrintWindow(h, dc, 2);
            g.ReleaseHdc(dc);
            bmp.Save(path, ImageFormat.Png);
        }
    }
}

public static class Drive
{
    public static System.Threading.Thread PressAsync(IntPtr h, string want)
    {
        var t = new System.Threading.Thread(() =>
        {
            try { Msaa.PressName(h, want); } catch (Exception e) { Console.Error.WriteLine("press: " + e.Message); }
        });
        t.IsBackground = true;
        t.SetApartmentState(System.Threading.ApartmentState.STA);
        t.Start();
        return t;
    }

    public static IntPtr WaitNew(int pid, List<IntPtr> known, int ms)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        while (sw.ElapsedMilliseconds < ms)
        {
            foreach (var w in Msaa.TopWindows(pid))
                if (!known.Contains(w)) return w;
            System.Threading.Thread.Sleep(100);
        }
        return IntPtr.Zero;
    }

    public static IntPtr Child(IntPtr top, int id)
    {
        foreach (var c in Msaa.Children(top))
            if (Msaa.CtrlId(c) == id) return c;
        return IntPtr.Zero;
    }
}
