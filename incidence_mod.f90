module incidence_mod
  implicit none
  private
  public :: incidence_matrix, read_edgelist, build_incidence, free_incidence
  public :: incidence_times_transpose_dense, column_sums

  integer, parameter :: dp = kind(1.0d0)

  !> 境界行列(接続行列) B (V x E) の疎行列表現 (CSC: 列圧縮格納)
  !> k列目 = k番目の辺。各列の非零要素は ちょうど2個 (始点: -1, 終点: +1)
  type :: incidence_matrix
    integer :: nv = 0                       ! 頂点数 V (行数)
    integer :: ne = 0                       ! 辺数   E (列数)
    integer, allocatable :: col_ptr(:)      ! 列kの要素は col_ptr(k) .. col_ptr(k+1)-1
    integer, allocatable :: row_idx(:)      ! 各非零要素の行番号(頂点番号)
    integer, allocatable :: val(:)          ! 各非零要素の値 (+1 or -1)
  end type incidence_matrix

contains

  !> エッジリストファイルを読み込む (空白/タブ区切り, CRLF可)
  subroutine read_edgelist(filename, edges, nv)
    character(len=*), intent(in)  :: filename
    integer, allocatable, intent(out) :: edges(:,:)   ! edges(1,k)=始点, edges(2,k)=終点
    integer, intent(out) :: nv                         ! 最大頂点番号
    integer :: u, ios, ne, k, a, b

    open(newunit=u, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'read_edgelist: cannot open file'

    ! 1回目: 辺数を数える
    ne = 0
    do
      read(u, *, iostat=ios) a, b
      if (ios /= 0) exit
      ne = ne + 1
    end do

    allocate(edges(2, ne))
    rewind(u)
    nv = 0
    do k = 1, ne
      read(u, *) edges(1,k), edges(2,k)
      nv = max(nv, edges(1,k), edges(2,k))
    end do
    close(u)
  end subroutine read_edgelist

  !> エッジリストから境界行列Bを構築 (k番目の辺 (u,v): 行u に -1, 行v に +1)
  subroutine build_incidence(edges, nv, B)
    integer, intent(in) :: edges(:,:)
    integer, intent(in) :: nv
    type(incidence_matrix), intent(out) :: B
    integer :: k, ne

    ne = size(edges, 2)
    B%nv = nv
    B%ne = ne
    allocate(B%col_ptr(ne+1), B%row_idx(2*ne), B%val(2*ne))

    do k = 1, ne
      if (edges(1,k) == edges(2,k)) stop 'build_incidence: self-loops are not supported'
      B%col_ptr(k) = 2*k - 1
      B%row_idx(2*k-1) = edges(1,k);  B%val(2*k-1) = -1
      B%row_idx(2*k)   = edges(2,k);  B%val(2*k)   = +1
    end do
    B%col_ptr(ne+1) = 2*ne + 1
  end subroutine build_incidence

  !> 各列の和 (正しい境界行列なら全て0)
  function column_sums(B) result(s)
    type(incidence_matrix), intent(in) :: B
    integer :: s(B%ne)
    integer :: k
    do k = 1, B%ne
      s(k) = sum(B%val(B%col_ptr(k):B%col_ptr(k+1)-1))
    end do
  end function column_sums

  !> 検証用: B * B^T (V x V, = グラフラプラシアン L = D - A) を密行列で計算
  function incidence_times_transpose_dense(B) result(L)
    type(incidence_matrix), intent(in) :: B
    integer :: L(B%nv, B%nv)
    integer :: k, p, q

    L = 0
    do k = 1, B%ne
      do p = B%col_ptr(k), B%col_ptr(k+1)-1
        do q = B%col_ptr(k), B%col_ptr(k+1)-1
          L(B%row_idx(p), B%row_idx(q)) = L(B%row_idx(p), B%row_idx(q)) &
                                          + B%val(p) * B%val(q)
        end do
      end do
    end do
  end function incidence_times_transpose_dense

  subroutine free_incidence(B)
    type(incidence_matrix), intent(inout) :: B
    if (allocated(B%col_ptr)) deallocate(B%col_ptr)
    if (allocated(B%row_idx)) deallocate(B%row_idx)
    if (allocated(B%val))     deallocate(B%val)
    B%nv = 0
    B%ne = 0
  end subroutine free_incidence

end module incidence_mod