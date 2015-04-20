Module Util
  Use set_precision
Contains
  Function cputime()
    Integer(wp) :: count, count_rate
    Real (wp) :: cputime
    Call System_clock(count,count_rate)
    cputime = Real(count,wp)/Real(count_rate,wp)
  End Function cputime
End Module Util
