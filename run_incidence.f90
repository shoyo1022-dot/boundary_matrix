program run_incidence
  use incidence_mod
  implicit none

  character(len=256) :: fname, outname
  integer, allocatable :: edges(:,:)
  integer, allocatable :: L(:,:)
  integer, allocatable :: Bd(:,:)
  type(incidence_matrix) :: B
  integer :: nv, k, i, j, ndiff, p, uo, ios, nargs
  integer, allocatable :: deg(:)

  ! --- ファイル名の取得 (引数 -> なければ対話入力) ---
  nargs = command_argument_count()

  if (nargs >= 1) then
    call get_command_argument(1, fname)
  else
    write(*, '(a)', advance='no') 'Input edge list file name: '
    read(*, '(a)') fname
  end if

  if (nargs >= 2) then
    call get_command_argument(2, outname)
  else
    write(*, '(a)', advance='no') 'Output file name (blank = incidence_out.dat): '
    read(*, '(a)') outname
  end if

  fname   = adjustl(fname)
  outname = adjustl(outname)
  if (len_trim(fname) == 0) stop 'run_incidence: no input file name given'
  if (len_trim(outname) == 0) outname = 'incidence_out.dat'

  call read_edgelist(trim(fname), edges, nv)
  call build_incidence(edges, nv, B)

  print '(a,i0,a,i0)', 'V = ', B%nv, ',  E = ', B%ne

  ! --- 疎行列(CSC)から V x E の密行列 Bd を作る ---
  allocate(Bd(B%nv, B%ne))
  Bd = 0
  do k = 1, B%ne
    do p = B%col_ptr(k), B%col_ptr(k+1)-1
      Bd(B%row_idx(p), k) = B%val(p)
    end do
  end do

  ! --- 出力ファイルへ V x E 行列を書き出す (1行 = 頂点1つ分) ---
  open(newunit=uo, file=trim(outname), status='replace', action='write', iostat=ios)
  if (ios /= 0) stop 'run_incidence: cannot open output file'

  do i = 1, B%nv
    write(uo, '(*(1x,i2))') (Bd(i,k), k = 1, B%ne)
  end do
  close(uo)

  print '(a,i0,a,i0,a,a)', 'Boundary matrix B (', B%nv, ' x ', B%ne, ') written to: ', trim(outname)

  ! --- 疎行列(CSC)の中身を表示 (先頭の数列) ---
  print '(a)', '--- Columns of B (k: source->target | row:value) ---'
  do k = 1, min(B%ne, 5)
    print '(i3,a,i0,a,i0,a)', k, ': ', edges(1,k), ' -> ', edges(2,k), '  |'
    do p = B%col_ptr(k), B%col_ptr(k+1)-1
      print '(a,i0,a,sp,i0)', '      row ', B%row_idx(p), ' : ', B%val(p)
    end do
  end do

  ! --- 検証1: 各列に非零がちょうど2個, 和 = 0 ---
  if (all(count(Bd /= 0, dim=1) == 2)) then
    print '(a)', '[OK] Every column has exactly two non-zero entries'
  else
    print '(a)', '[NG] Some columns do not have exactly two non-zero entries'
  end if
  if (all(sum(Bd, dim=1) == 0)) then
    print '(a)', '[OK] All column sums are 0 (one +1 and one -1 per edge)'
  else
    print '(a)', '[NG] Some columns have a non-zero sum'
  end if

  ! --- 検証2: B B^T = D - A (ラプラシアン) ---
  allocate(L(nv,nv), deg(nv))
  L = incidence_times_transpose_dense(B)

  deg = 0
  do k = 1, B%ne
    deg(edges(1,k)) = deg(edges(1,k)) + 1
    deg(edges(2,k)) = deg(edges(2,k)) + 1
  end do

  ndiff = 0
  do i = 1, nv
    if (L(i,i) /= deg(i)) ndiff = ndiff + 1
    do j = 1, nv
      if (i == j) cycle
      if (L(i,j) > 0) ndiff = ndiff + 1          ! 非対角は 0 か 負
    end do
  end do
  if (ndiff == 0) then
    print '(a)', '[OK] Diagonal of B*B^T = degree, off-diagonal <= 0 (Laplacian)'
  else
    print '(a,i0)', '[NG] Number of mismatches: ', ndiff
  end if

  print '(a,i0)', 'Degree of vertex 1 (L(1,1) of B*B^T) = ', L(1,1)

  call free_incidence(B)
end program run_incidence