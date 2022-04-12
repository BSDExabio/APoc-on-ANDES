**************************************************************************
*     APoc is a program for assessing protein pocket structural        *
*     similarity. The source code evolves from iAlign and TM-align.      *
*                                                                        *
*     Mu Gao  <mu.gao@gatech.edu>                                        *
**************************************************************************
      module Constants
        implicit none
        include 'pars.h'
      end module

      module PDBType
        use Constants
        implicit none
        type :: pdbstruct
          !!! regular pdb structure
          character*500 :: filename
          character*3,pointer :: resname(:)
          character,  pointer :: seq(:),icode(:),chain(:)
          integer,pointer :: resid(:)
          integer,pointer :: sec_struct(:),iseq(:)
          real,pointer :: coor_ca(:,:),coor_cb(:,:)
          real,pointer :: cbvec(:,:)
          integer :: num_residue

          !!! pocket structure
          integer :: num_pocket
          integer,pointer :: pk_res(:,:),pk_len(:),pk_vol(:)  !pocket residues, length, volume
          character*30,pointer :: pk_name(:)
        end type pdbstruct
      end module

      module Options
        use Constants
        implicit none
        type :: alnopt
          integer score_flag    ! 0 - TM-score, 1 - IS-score ( contact-weighted TM-score ), 2 - PS-score
          integer order_flag    ! 0 - sequential order independent, 1 - dependent
          integer verbo_flag    ! 2 - print alignment in pairs, default 1 in sequential, 2 in non-sequential mode
          integer quick_flag    ! 0 - slow score caculations, 1 - regular score calculations (default)
          integer seqd0_flag    ! 0 - use non-sequential d0, 1 - use sequential d0 (for testing purpose)
          integer ds_flag       ! 0 - default, 1 - d0 search based on first struct, 2 - d0 search based on second struct
          integer skip_aln_flag ! 0 - default, 1 - skip alignment search, take the input pdb as aligned pairs, just calculate score
          integer prealn_flag   ! 0 - default, 1 - use pre-defined alignment as the starting seed to find the optimal alignment
          integer lt_file_flag  ! 0 - default, no list file defined for templates, 1 - read a list of templates from a file
          integer lq_file_flag  ! 0 - default, no list file defined for queries, 1 - read a list of queries from a file
          integer block_flag    ! 0 - default, no pdb block file to read, 1 - read a block pdb file that contains multiple pdb records
                                ! in one file, this could improve I/O dramatically

          !!! pocket related options
          real    pk_min_vol  ! minimal volume of pocket
          integer pk_min_len  ! minimal number of pocket residues
          integer npksel1,npksel2
          character*30 pksel1(maxp),pksel2(maxp) !selected pockets for comparison


          !!! old TM-align options
          integer m_out      !decided output
          integer m_fix      !fixed length-scale only for output
          integer m_ave      !using average length
          integer m_d0_min   !miminum d0 for search
          integer m_d0       !given d0 for both search and output
          integer L_fix

          character*500 outname,lt_file,lq_file,block_file
          real d0_min_input
          real d0_fix
        end type alnopt
      end module

      module AlnResult
        use Constants
        implicit none
        type :: alignment
          integer length
          integer map1(maxr),map2(maxr)
        end type

        type :: alnres
          real score
          real rmsd
          real seqid
          real d0
          real*8 pvalue,zscore

          integer seqsc,seqsc_pos
          integer init

          character*2 measure
          type(alignment) aln
        end type
      end module


**************************************************************************
**************************************************************************
*   Main Program
*
**************************************************************************
**************************************************************************
      program palign
      use Constants
      use PDBType
      use Options
      use AlnResult
      implicit none

      type(pdbstruct) :: struct1,struct2,struct(maxt)
      type(alnopt) :: opt,opt_tmp
      type(alnres) :: fullout,pkout
      type(alignment) :: prealn

      character*500 fnam,pdb(100),templ_lst(maxt),query_lst(maxt)
      character*20 numstr
      character*30 pname
      character*30 lt_pksel(maxt,maxp),lq_pksel(maxt,maxp) !selected pockets for comparison
      integer i,j,pkind1,pkind2,narg,iargc
      integer npk1,npk2,pk1(maxp),pk2(maxp)
      integer lt_npksel(maxt), lq_npksel(maxt)
      integer templ_str_ind(maxt),query_str_ind(maxt)
      integer templ_num,query_num,struct_num,itempl,iquery
      integer tot_num_pair
      integer is_pk_found

******* options ----------->
      opt%m_out=-1                  !decided output
      opt%m_fix=-1                  !fixed length-scale only for output
      opt%m_ave=-1                  !using average length
      opt%m_d0_min=-1               !diminum d0 for search
      opt%m_d0=-1                   !given d0 for both search and output
      opt%score_flag=2        !0 - TM-score, 1- IS-score, 2 - PS-score (default)
      opt%order_flag=0        !0 - order independent, 1- order dependent
      opt%seqd0_flag=0        !0 - normal, 1 - use sequential d0 in non-sequence aln (for testing purpose)
      opt%verbo_flag=1
      opt%quick_flag=1
      opt%skip_aln_flag=0
      opt%prealn_flag=1
      opt%ds_flag=0
      opt%lt_file_flag=0
      opt%lq_file_flag=0
      opt%block_flag=0

      opt%npksel1=0
      opt%npksel2=0
      opt%pk_min_vol=100
      opt%pk_min_len=10

      templ_num=0
      query_num=0
      tot_num_pair=0

      narg=iargc()
      i=0
      j=0
 115  continue
      i=i+1
      call getarg(i,fnam)
      if(fnam.eq.' '.or.fnam.eq.'?'.or.fnam.eq.'-h')then
         call printHelp()
         goto 9999
      else if(fnam.eq.'-L')then  !change both L_all and d0
         opt%m_fix=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)opt%L_fix
      else if(fnam.eq.'-dmin')then
         opt%m_d0_min=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)opt%d0_min_input
      else if(fnam.eq.'-d0')then
         opt%m_d0=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)opt%d0_fix
      else if(fnam.eq.'-a')then ! this will change the superposed output but not the alignment
        opt%m_ave=1
      elseif(fnam.eq.'-b')then
        opt%m_ave=2
      else if(fnam.eq.'-c')then
        opt%m_ave=3
      else if(fnam.eq.'-m')then
         i=i+1
         call getarg(i,fnam)
         if(fnam.eq.'tm')then
            opt%score_flag=0
         elseif(fnam.eq.'is')then
            opt%score_flag=1
         elseif(fnam.eq.'ps')then
            opt%score_flag=2
         else
            print '(3A)','Error: ',fnam, ' is not a valid option'
            call printHelp()
            goto 9999
         endif
      else if(fnam.eq.'-sod')then
         opt%order_flag=1
      else if(fnam.eq.'-skipaln')then
         opt%skip_aln_flag=1
      else if(fnam.eq.'-seqd0')then ! use sequential d0 for testing purpose
         opt%seqd0_flag=1
      else if(fnam.eq.'-v')then ! verbose option
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)opt%verbo_flag
         if(opt%verbo_flag.lt.0.or.opt%verbo_flag.gt.2)then
            call printHelp()
            goto 9999
         endif
      else if(fnam.eq.'-q')then ! quick mode
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)opt%quick_flag
         if(opt%quick_flag.lt.0.or.opt%quick_flag.gt.1)then
            call printHelp()
            goto 9999
         endif
      else if(fnam.eq.'-ds')then !
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)opt%ds_flag
         if(opt%ds_flag.lt.0.or.opt%ds_flag.gt.2)then
            call printHelp()
            goto 9999
         endif
      elseif(fnam.eq.'-fa')then
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)opt%prealn_flag
         if(opt%prealn_flag.ne.0.and.opt%prealn_flag.ne.1)then
            call printHelp()
            goto 9999
         endif
      elseif(fnam.eq.'-pvol')then
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)opt%pk_min_vol
      elseif(fnam.eq.'-plen')then
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)opt%pk_min_len
      elseif(fnam.eq.'-pt')then
         i=i+1
         call getarg(i,fnam)
         call getPKind(fnam,opt%pksel1,opt%npksel1)
      elseif(fnam.eq.'-pq')then
         i=i+1
         call getarg(i,fnam)
         call getPKind(fnam,opt%pksel2,opt%npksel2)
      elseif(fnam.eq.'-lt')then
         i=i+1
         call getarg(i,fnam)
         read(fnam,'(A)')opt%lt_file
         opt%lt_file_flag=1
         call readListFile(opt%lt_file,templ_lst,templ_num,lt_npksel,lt_pksel)
         if(templ_num.eq.0)then
            print *,'Error: no structure defined in file: ',opt%lt_file
            call printHelp()
            goto 9999
         endif
      elseif(fnam.eq.'-lq')then
         i=i+1
         call getarg(i,fnam)
         read(fnam,'(A)')opt%lq_file
         opt%lq_file_flag=1
         call readListFile(opt%lq_file,query_lst,query_num,lq_npksel,lq_pksel)
         if(templ_num.eq.0)then
            print *,'Error: no structure defined in file: ',opt%lq_file
            call printHelp()
            goto 9999
         endif
      elseif(fnam.eq.'-block')then
         i=i+1
         call getarg(i,fnam)
         read(fnam,'(A)')opt%block_file
         opt%block_flag=1
      elseif(fnam(1:1).eq.'-')then
         print *,'Error: Invalid option ',trim(adjustl(fnam))
         call printHelp()
         goto 9999
      else
         j=j+1
         if(j.gt.2)then
            print *, 'Error: unknown parameter ',fnam
            call printHelp()
            goto 9999
         endif
         pdb(j)=fnam
      endif
      if(i.lt.narg)goto 115


      write(*,*)
      write(*,*)'*****************************************************',
     &     '*********************'
      write(*,*)'                                APoc                 '

      write(*,*)'             Structural alignment of protein pockets ',
     &     '                 '
      write(*,*)'*****************************************************',
     &     '*********************'
      write(*,*)


      if(opt%lt_file_flag.eq.0.and.opt%lq_file_flag.eq.0)then
         templ_lst(1)=pdb(1)
         query_lst(1)=pdb(2)

         templ_num=1
         query_num=1
      else if(opt%lq_file_flag.eq.0)then
         do i=1,j
            query_lst(i)=pdb(j)
         enddo
         query_num=j
      endif


      struct_num=0
      !!! read block file first if specified
      print '(A)', 'Loading PDB files, this could take a while ...'
      if(opt%block_flag.eq.1)then
         call readPDBBlock( opt%block_file, struct_num, struct )
      endif

      !!! read the structures of templates
      do i=1,templ_num
         templ_str_ind(i)=0
         do j=1,struct_num
            if(templ_lst(i).eq.struct(j)%filename)then
               templ_str_ind(i)=j
               exit
            endif
         enddo
         if(templ_str_ind(i).eq.0)then
            if(struct_num.eq.maxt)then
               print *,'Error: exceed the maximum of structures allowed ',maxt
               goto 9999
            endif
            struct_num=struct_num+1
            call readPDBFileNew( templ_lst(i), struct(struct_num) )
            templ_str_ind(i)=struct_num
         endif
      enddo

      !!! read the structures of queries
      do i=1,query_num
         query_str_ind(i)=0
         do j=1,struct_num
            if(query_lst(i).eq.struct(j)%filename)then
               query_str_ind(i)=j
               exit
            endif
         enddo
         if(query_str_ind(i).eq.0)then
            if(struct_num.eq.maxt)then
               print *,'Error: exceed the maximum of structures allowed ',maxt
               goto 9999
            endif

            struct_num=struct_num+1
            call readPDBFileNew( query_lst(i), struct(struct_num) )
            query_str_ind(i)=struct_num
         endif
      enddo
      print '(A, /)', 'PDB files loaded'




      write(*,'(A,I7,A)')'Minimal volume of pocket: ',int(opt%pk_min_vol), ' grid pts'
      write(*,'(A,I7,A)')'Minimal length of pocket: ',opt%pk_min_len, ' residues'
      write(*,'(A,I5,A)')'Total number of template:   ',templ_num, ' structure(s)'
      write(*,'(A,I5,A)')'Total number of query:      ',query_num, ' structure(s)'
      write(*,*)





cccccc loop query/targets
      do 10 iquery=1,query_num
      struct2 = struct(query_str_ind(iquery))

cccccccccc Select qualifying pockets   ccccccccccccccccccccccc
      npk2=0
      if(opt%npksel2.eq.0.and.lq_npksel(iquery).eq.0)then
         opt%npksel2=struct2%num_pocket
         do i=1,struct2%num_pocket
            opt%pksel2(i)=struct2%pk_name(i)
         enddo
      else if(lq_npksel(iquery).gt.0)then
         opt%npksel2=0
         do j=1,lq_npksel(iquery)
            is_pk_found=0
            do i=1,struct2%num_pocket
               pname=lq_pksel(iquery,j)
               if(pname.eq.struct2%pk_name(i))then
                  opt%pksel2(j)=struct2%pk_name(i)
                  opt%npksel2=opt%npksel2+1
                  is_pk_found=is_pk_found+1
                  exit
               endif
            enddo

            if(is_pk_found.eq.0)then
               print '(4A)', 'Warning: pocket ',trim(pname),' not found in ',trim(query_lst(iquery))
c               goto 9999
            endif
         enddo
      endif

      do i=1,opt%npksel2
         pname=opt%pksel2(i)
         is_pk_found=0
         do j=1,struct2%num_pocket
            if(pname.eq.struct2%pk_name(j))then
               is_pk_found=is_pk_found+1
               exit
            endif
         enddo

         if(is_pk_found.eq.0)then
            print '(4A)', 'Warning: pocket ',trim(pname),' not found in ',trim(query_lst(iquery))
c            goto 9999
         endif

         if(struct2%pk_vol(j).ge.opt%pk_min_vol.and.
     &        struct2%pk_len(j).ge.opt%pk_min_len) then
            npk2=npk2+1
            pk2(npk2)=j
         endif
      enddo


cccccc loop templates
      do 20 itempl=1,templ_num
      struct1 = struct(templ_str_ind(itempl))

cccccccccc Select qualifying pockets   ccccccccccccccccccccccc
      npk1=0
      if(opt%npksel1.eq.0.and.lt_npksel(itempl).eq.0)then
         opt%npksel1=struct1%num_pocket
         do i=1,struct1%num_pocket
            opt%pksel1(i)=struct1%pk_name(i)
         enddo
      else if(lt_npksel(itempl).gt.0)then
         opt%npksel1=0
         do j=1,lt_npksel(itempl)
            is_pk_found=0
            do i=1,struct1%num_pocket
               pname=lt_pksel(itempl,j)
               if(pname.eq.struct1%pk_name(i))then
                  opt%pksel1(j)=struct1%pk_name(i)
                  opt%npksel1=opt%npksel1+1
                  is_pk_found=is_pk_found+1
                  exit
               endif
            enddo

            if(is_pk_found.eq.0)then
               print '(4A)', 'Warning: pocket ',trim(pname),' not found in ',trim(templ_lst(itempl))
c               goto 9999
            endif
         enddo
      endif

      do i=1,opt%npksel1
         pname=opt%pksel1(i)
         is_pk_found=0
         do j=1,struct1%num_pocket
            if(pname.eq.struct1%pk_name(j))then
               is_pk_found=is_pk_found+1
               exit
            endif
         enddo

         if(is_pk_found.eq.0)then
            print '(4A)', 'Warning: pocket ',trim(pname),' not found in ',trim(templ_lst(itempl))
            call printHelp()
c            goto 9999
         endif

         if(struct1%pk_vol(j).ge.opt%pk_min_vol.and.
     &        struct1%pk_len(j).ge.opt%pk_min_len) then
            npk1=npk1+1
            pk1(npk1)=j
         endif
      enddo

      tot_num_pair=tot_num_pair+1
      write(numstr,'(I)')tot_num_pair
      write(*,'(3A)')'============== Structure Pair ', trim(adjustl(numstr)), " =============="
      write(*,*)

      write(*,'(A,I3,2A)')'Found ',npk1,' qualifying pockets in ',trim(templ_lst(itempl))
      write(*,'(A,I3,2A)')'Found ',npk2,' qualifying pockets in ',trim(query_lst(iquery))
      write(*,*)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccc Run full chain alignment, this is a regular TM-alignment
      if(opt%prealn_flag.eq.1)then
         opt_tmp = opt
         opt_tmp%prealn_flag = 0
         opt_tmp%order_flag = 1
         opt_tmp%score_flag = 0
         prealn%length = 0
         write(*,'(A)')'>>>>>>>>>>>>>>>>>>>>>>>>>   '//
     &        'Global alignment   <<<<<<<<<<<<<<<<<<<<<<<<<<'
         call comp2str( struct1,0,struct2,0,opt_tmp,fullout,prealn )
         write(*,*)
      endif

ccccccc Run pocket alignment
      opt_tmp = opt
      opt_tmp%prealn_flag = 1
      do j = 1,npk2
         pkind2 = pk2(j)
         do i = 1,npk1
            pkind1 = pk1(i)
            write(*,'(A)')'>>>>>>>>>>>>>>>>>>>>>>>>>   '//
     &       'Pocket alignment   <<<<<<<<<<<<<<<<<<<<<<<<<<'
            call mkPreAln( struct1,pkind1,struct2,pkind2,fullout,prealn )
            call comp2str( struct1,pkind1,struct2,pkind2,opt_tmp,pkout,prealn)
            write(*,*)
         enddo
      end do

 20   continue
 10   continue

 9999 do i=1,struct_num
         deallocate(struct(i)%resname)
         deallocate(struct(i)%seq)
         deallocate(struct(i)%icode)
         deallocate(struct(i)%chain)
         deallocate(struct(i)%resid)
         deallocate(struct(i)%sec_struct)
         deallocate(struct(i)%iseq)
         deallocate(struct(i)%coor_ca)
         deallocate(struct(i)%coor_cb)
         deallocate(struct(i)%cbvec)
         deallocate(struct(i)%pk_res)
         deallocate(struct(i)%pk_len)
         deallocate(struct(i)%pk_vol)
         deallocate(struct(i)%pk_name)
      enddo

      END






***********************************************************************
***********************************************************************
*     Comparison of two structures
***********************************************************************
***********************************************************************
      subroutine comp2str(struct1,pkind1,struct2,pkind2,opt,out,prealn)
      use Constants
      use PDBType
      use Options
      use AlnResult

      type(pdbstruct) struct1,struct2
      type(alnopt) :: opt
      type(alnres) :: out
      type(alignment) :: prealn

      integer pkind1, pkind2

      COMMON/BACKBONE/XA(3,maxr,0:1)
      common/dpc/score(maxr,maxr),gap_open

      common/cbvec/cbvec(3,maxr,0:1)  ! normalized vector from ca to cb
      common/cbeta/xtb(maxr),ytb(maxr),ztb(maxr)

      integer invmap0(maxr)     !final alignment from second sequence to first
      integer invmap0_r(maxr)   !final alignment from first  sequence to second

      common/length/nseq1,nseq2
      common/d0/d0,anseq
      common/d0min/d0_min
      common/d00/d00,d002

      character*500 fnam,outname,file1,file2
      character*500 pdbfile1,pdbfile2
      character*3 aa(-1:20),aanam,ss1(maxr),ss2(maxr)
      character*100 s,du
      character*200 outnameall_tmp,outnameall
      character seq1(0:maxr),seq2(0:maxr)
      character chain1(maxr),chain2(maxr)
      character icode1(maxr),icode2(maxr)
      character aseq1(maxr2),aseq2(maxr2),aseq3(maxr2)
      character*2   measure   ! 'TM' -> TM-score, 'IS' -> IS-score
      character*3  note
      character*20 numstr,numstr2

      dimension xx(maxr),yy(maxr),zz(maxr)
      dimension m1(maxr),m2(maxr)
      dimension xtm1(maxr),ytm1(maxr),ztm1(maxr)
      dimension xtm2(maxr),ytm2(maxr),ztm2(maxr)
      common/init/invmap_i(maxr)

      common/TM/TM,TMmax
      common/n1n2/n1(maxr),n2(maxr)
      common/d8/d8
      common/initial4/mm1(maxr),mm2(maxr)

      real  contcf   ! Ca-Ca contact cutoff
      integer cont1(maxc,maxr),cont2(maxc,maxr)     !contact list
      common/contlst/cont1,cont2
      common/contnum/ncont1(maxr),ncont2(maxr)      !number of contacts
      common/contdis/dcont1,dcont2                  !contact distance
      real  dcont1(maxr,maxr),dcont2(maxr,maxr)
      common/contcf/d_col,d_col2

      common/sec/isec(maxr),jsec(maxr)

      integer is_ali(maxr, maxr)   !alignment map
      integer col,ncol
      integer best_mode,best_init

      real    fcol
      real    op   !order paramter
      real    occupancy

ccc   RMSD:
      double precision r_1(3,maxr),r_2(3,maxr),r_3(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /maxr*1.0/
ccc

      integer seqsc,tot_seqsc,tot_seqsc_pos
      integer seqmat(0:19,0:19), blmat(0:19,0:19)    ! sequence substitution matrix, e.g., BLOSUM62
      integer iseq1(maxr),iseq2(maxr)
      real    sdel   ! penalty scaling factor for unmatched amino acids
      common/seq/iseq1,iseq2,seqmat,sdel
      common/blmat/blmat

      real*8   pvalue, best_pvalue, zs, best_zs
      real     eps

ccc   Timing:
      real*8 total_time
      integer ( kind = 4 ) clock_start
      integer ( kind = 4 ) clock_stop
      integer ( kind = 4 ) clock_max
      integer ( kind = 4 ) clock_rate



      call system_clock ( clock_start, clock_rate, clock_max )


ccccc Initialization
      pdbfile1 = struct1%filename
      pdbfile2 = struct2%filename

      nseq1 = struct1%num_residue
      if(pkind1.gt.0) nseq1 = struct1%pk_len(pkind1)
      do i = 1, nseq1
         k = i
         if(pkind1.gt.0) k = struct1%pk_res(pkind1,i)

         chain1(i) = struct1%chain(k)
         icode1(i) = struct1%icode(k)
         iseq1(i) = struct1%iseq(k)
         isec(i) = struct1%sec_struct(k)
         seq1(i) = struct1%seq(k)
         mm1(i) = struct1%resid(k)
         ss1(i) = struct1%resname(k)

         do j=1,3
            xa(j,i,0)=struct1%coor_ca(j,k)
            cbvec(j,i,0)=struct1%cbvec(j,k)
         enddo
      enddo

      nseq2 = struct2%num_residue
      if(pkind2.gt.0) nseq2 = struct2%pk_len(pkind2)
      do i = 1, nseq2
         k = i
         if(pkind2.gt.0) k = struct2%pk_res(pkind2,i)

         chain2(i) = struct2%chain(k)
         icode2(i) = struct2%icode(k)
         iseq2(i) = struct2%iseq(k)
         jsec(i) = struct2%sec_struct(k)
         seq2(i) = struct2%seq(k)
         mm2(i) = struct2%resid(k)
         ss2(i) = struct2%resname(k)

         do j=1,3
            xa(j,i,1)=struct2%coor_ca(j,k)
            cbvec(j,i,1)=struct2%cbvec(j,k)
         enddo
      enddo
ccccc End of initialization



c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      contcf=15
      call calCaCont(xa,0,nseq1,ncont1,cont1,contcf,dcont1)  !contacts of structure 1
      call calCaCont(xa,1,nseq2,ncont2,cont2,contcf,dcont2)  !contacts of structure 2

      !do i=1,nseq1
      !   print *,i,ncont1(i),(cont1(j,i),j=1,ncont1(i))
      !enddo


      measure='TM'
      if(opt%order_flag.eq.0.and.opt%verbo_flag.eq.1)opt%verbo_flag=2
      if(opt%score_flag.eq.1)then
         measure='IS'
      else if(opt%score_flag.eq.2)then
         measure='PS'
      endif

      sdel=0.8
      call readAA8(seqmat)
      !call readUniMat(seqmat)
      call readBlosum62(blmat)


*!!!  Define parameters for searching the optimal alignment --------->
      d0_min=1.0
      if(opt%m_d0_min.eq.1)then
         d0_min=opt%d0_min_input    !for search
      endif
      anseq_min=min(nseq1,nseq2)
      anseq=anseq_min              !length for defining TMscore in search
      if(opt%ds_flag.eq.1)anseq=nseq1
      if(opt%ds_flag.eq.2)anseq=nseq2

      d8=1.5*anseq_min**0.3+3.5 !remove pairs with dis>d8 during search & final
      if(opt%order_flag.eq.0.and.d8.gt.8) d8=8

      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif

      if((opt%order_flag.eq.0.and.opt%seqd0_flag.eq.0).or.opt%score_flag.eq.2)then
         if(anseq.gt.5)d0=0.7*(anseq-5)**(0.25)-0.2
         if(d0.lt.0.5)d0=0.5
      endif


      if(d0.lt.d0_min)d0=d0_min
      if(opt%m_d0.eq.1)d0=opt%d0_fix
      d00=d0                    !for quickly calculate TM-score in searching
      if(d00.gt.8)d00=8
      if(d00.lt.4.5)d00=4.5
      d002=d00**2

      d_col=d8          !distance cutoff for aligned residues with contact overlap
      d_col2=d_col**2   !prevent costs of sqrt


      nseq=max(nseq1,nseq2)
      do i=1,nseq
         n1(i)=i
         n2(i)=i
      enddo

      zs=-1


***** find the optimal alignment **************************
      if(opt%skip_aln_flag.eq.0)then
         CALL super_align(invmap0,opt%score_flag,opt%order_flag,
     &      opt%quick_flag,best_init,best_mode,prealn)                 !to find invmap from seq2 to seq1
      else
         do j=1,nseq2
            if(j.le.nseq1)then
               invmap0(j)=j     ! no alignment search, structures are aligned
            else
               invmap0(j)=-1
            endif
         enddo
      endif


************************************************************
***   resuperpose to find residues of dis<d8 ------------------------>
      n_al=0
      if(best_mode.eq.0)then
         do j=1,nseq2
            if(invmap0(j).gt.0)then
               i=invmap0(j)
               n_al=n_al+1
               xtm1(n_al)=xa(1,i,0)
               ytm1(n_al)=xa(2,i,0)
               ztm1(n_al)=xa(3,i,0)
               xtm2(n_al)=xa(1,j,1)
               ytm2(n_al)=xa(2,j,1)
               ztm2(n_al)=xa(3,j,1)
               m1(n_al)=i       !for recording residue order
               m2(n_al)=j
            endif
         enddo
      else
         do i=1,nseq1
            invmap0_r(i)=-1     !seq1 to seq2, reverse of invmap
         enddo

         do j=1,nseq2
            i=invmap0(j)
            if(i.gt.0) invmap0_r(i)=j
         enddo

         do i=1,nseq1
            j=invmap0_r(i)      !j aligned to i
            if(j.gt.0)then
               n_al=n_al+1
               xtm1(n_al)=xa(1,i,0) !for TM-score
               ytm1(n_al)=xa(2,i,0)
               ztm1(n_al)=xa(3,i,0)
               xtm2(n_al)=xa(1,j,1)
               ytm2(n_al)=xa(2,j,1)
               ztm2(n_al)=xa(3,j,1)

               m1(n_al)=i       !aligned position to original index
               m2(n_al)=j
            endif
         enddo
      endif


      d0_input=d0
      d0_search=d0
      if (d0_search.gt.8.0)d0_search=8.0
      if (d0_search.lt.1.0)d0_search=1.0
      call TMsearch(d0_input,d0_search,n_al,xtm1,ytm1,ztm1,m1,n_al,
     &     xtm2,ytm2,ztm2,m2,TM,Rcomm,Lcomm,2,opt%score_flag,ncol) !calcualte TMscore with dis<d8 only


***   remove dis>d8 in normal TM-score calculation for final report----->
***   also remove aligned residues with zero contact overlap for IS-score
      call init_ali_map( is_ali, nseq1, nseq2 )
      do k=1,n_al
         i=m1(k)
         j=m2(k)
         is_ali(i,j) = 1
      enddo

      j=0
      n_eq=0
      tot_seqsc=0
      tot_seqsc_pos=0
      do i=1,n_al
         dis2=sqrt((xtm1(i)-xtm2(i))**2+(ytm1(i)-ytm2(i))**2+
     &        (ztm1(i)-ztm2(i))**2)

         if(dis2.le.d8)then
            j=j+1
            xtm1(j)=xa(1,m1(i),0) !for TM-score
            ytm1(j)=xa(2,m1(i),0)
            ztm1(j)=xa(3,m1(i),0)
            xtm2(j)=xa(1,m2(i),1)
            ytm2(j)=xa(2,m2(i),1)
            ztm2(j)=xa(3,m2(i),1)

            m1(j)=m1(i)
            m2(j)=m2(i)
            if(ss1(m1(i)).eq.ss2(m2(i)))then
               n_eq=n_eq+1
            endif
            iaa1=iseq1(m1(i))
            iaa2=iseq2(m2(i))
            if(iaa1.ge.0.and.iaa2.ge.0)then
               seqsc=blmat(iaa1,iaa2)
               tot_seqsc=tot_seqsc+seqsc
               if(seqsc.gt.0)tot_seqsc_pos=tot_seqsc_pos+seqsc
            endif
         endif
      enddo

      n8_al=j
      seq_id=real(n_eq)/real(n8_al) !sequence identity over aligned region

ccccc !for calculating p-value of scores, always use short sequence length to derive d0_search
ccccc !recalculate the score is necessary since aligned residues are slightly different by excluding d > d8
      d0_input=d0
      d0_search=d0
      if (d0_search.gt.8.0)d0_search=8.0
      if (d0_search.lt.1.0)d0_search=1.0
      call TMsearch(d0_input,d0_search,n8_al,xtm1,ytm1,ztm1,m1,n8_al,
     &     xtm2,ytm2,ztm2,m2,TM8,Rcomm,Lcomm,3,opt%score_flag,ncol) !normal TMscore

      TM8=TM8*n8_al/anseq

      pvalue=1
      if( opt%score_flag.eq.2 )then
         f0 = 0.23 - 12*(anseq**(-1.88))
         TM8=(TM8+f0)/(1+f0)      ! rescale PS-score
         call calcPvalue(TM8,nseq1,nseq2,pvalue,zs)
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




*!!!  Final score output is based on the second protein------------------>
      d0_min=0.5                !for output
      anseq=nseq2               !length for defining final TMscore
      if(opt%m_ave.eq.1)anseq=(nseq1+nseq2)/2.0 !<L>
      if(opt%m_ave.eq.2)anseq=min(nseq1,nseq2)
      if(opt%m_ave.eq.3)anseq=max(nseq1,nseq2)
      if(anseq.lt.anseq_min)anseq=anseq_min
      if(opt%m_fix.eq.1)anseq=opt%L_fix !input length

      dslen=anseq
      if(opt%ds_flag.eq.1)dslen=nseq1
      if(opt%ds_flag.eq.2)dslen=nseq2

      if(dslen.gt.15)then
         d0=1.24*(dslen-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif

      if((opt%order_flag.eq.0.and.opt%seqd0_flag.eq.0).or.opt%score_flag.eq.2)then
         if(dslen.gt.5)d0=0.7*(dslen-5)**(0.25)-0.2
         if(d0.lt.0.5)d0=0.5
      endif

      if(d0.lt.d0_min)d0=d0_min
      if(opt%m_d0.eq.1)d0=opt%d0_fix

      !final score
      d0_input=d0
      d0_search=d0
      if (d0_search.gt.8.0)d0_search=8.0
      if (d0_search.lt.1.0)d0_search=1.0

      do i=1,n8_al
         xtm1(i)=xa(1,m1(i),0)  !for PS-score, important for rotating Ca-Cb vectors properly
         ytm1(i)=xa(2,m1(i),0)
         ztm1(i)=xa(3,m1(i),0)
         xtm2(i)=xa(1,m2(i),1)
         ytm2(i)=xa(2,m2(i),1)
         ztm2(i)=xa(3,m2(i),1)
      enddo

      call TMsearch(d0_input,d0_search,n8_al,xtm1,ytm1,ztm1,m1,n8_al,
     &     xtm2,ytm2,ztm2,m2,TM8,Rcomm,Lcomm,3,opt%score_flag,ncol) !normal TMscore
      rmsd8_al=Rcomm
      TM8=TM8*n8_al/anseq       !TM-score after cutoff

      if(opt%score_flag.eq.2)then
         f0 = 0.23 - 12*(anseq**(-1.88))
         TM8=(TM8+f0)/(1+f0)      ! rescale PS-score
      endif


      out%score = TM8
      out%measure = measure
      out%seqid = seq_id
      out%rmsd  = rmsd8_al
      out%pvalue= pvalue
      out%zscore= zs
      out%d0 = d0
      out%init = best_init
      out%seqsc = tot_seqsc
      out%seqsc_pos = tot_seqsc_pos

      out%aln%length = n8_al
      do i=1,n8_al
         out%aln%map1(i) = m1(i)
         out%aln%map2(i) = m2(i)
      enddo



********* for output summary ******************************
      call getBasename( pdbfile1, file1 )
      call getBasename( pdbfile2, file2 )

      if( pkind1.eq.0 )then
         write(*,101)1,file1,nseq1
      else
         write(*,102)1,file1,nseq1,trim(struct1%pk_name(pkind1))
      endif
      if( pkind2.eq.0 )then
         write(*,101)2,file2,nseq2
      else
         write(*,102)2,file2,nseq2,trim(struct2%pk_name(pkind2))
      endif
 101  format('Structure ',I1,': ',A30,'  Length =',I4,' AAs,', ' Full')
 102  format('Structure ',I1,': ',A30,'  Length =',I4,' AAs,', ' Pocket:',A)
      write(*,*)

      write(numstr, '(I)')n8_al
      if(opt%score_flag.eq.2)then
         write(*,103)measure,TM8,pvalue,zs
         write(*,155)adjustl(numstr),rmsd8_al,seq_id
      else
         write(*,154)measure,TM8
         write(*,155)adjustl(numstr),rmsd8_al,seq_id
      endif
 103  format(A2,'-score =',F8.5,', P-value =',E12.4E3,', Z-score =',F7.3)
 154  format(A2,'-score =',F8.5)
 155  format('Number of aligned residues  = ', A, /,
     &       'RMSD =',F6.2, ', Seq identity  =',F6.3, /)

********* extract rotation matrix ------------>
      L=0
      do i=1,n8_al
         k=m1(i)
         L=L+1
         r_1(1,L)=xa(1,k,0)
         r_1(2,L)=xa(2,k,0)
         r_1(3,L)=xa(3,k,0)
         r_2(1,L)=xtm1(i)
         r_2(2,L)=ytm1(i)
         r_2(3,L)=ztm1(i)
       enddo
       if(L.ge.3)then
         call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2
         armsd=dsqrt(rms/L)
         write(*,*)'-------- rotation matrix to rotate Chain-1 to ',
     &        'Chain-2 ------'
         write(*,*)'i          t(i)         u(i,1)         u(i,2) ',
     &        '        u(i,3)'
         do i=1,3
            write(*,204)i,t(i),u(i,1),u(i,2),u(i,3)
         enddo
c         do i=1,nseq1
c            ax=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
c            ay=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
c            az=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
c         enddo
         write(*,*)
      endif
 204  format(I2,f18.10,f15.10,f15.10,f15.10)



      if(opt%verbo_flag.eq.1)then
************  output aligned sequences **************************
         ii=0
         i1_old=1
         i2_old=1
         do i=1,n8_al
            do j=i1_old,m1(i)-1
               ii=ii+1
               aseq1(ii)=seq1(j)
               aseq2(ii)='-'
               aseq3(ii)=' '
            enddo
            do j=i2_old,m2(i)-1
               ii=ii+1
               aseq1(ii)='-'
               aseq2(ii)=seq2(j)
               aseq3(ii)=' '
            enddo
            ii=ii+1
            aseq1(ii)=seq1(m1(i))
            aseq2(ii)=seq2(m2(i))
            dis2=sqrt((xtm1(i)-xtm2(i))**2+
     &           (ytm1(i)-ytm2(i))**2+(ztm1(i)-ztm2(i))**2)
            if(dis2.le.5)then
               aseq3(ii)=':'
            else
               aseq3(ii)='.'
            endif
            i1_old=m1(i)+1
            i2_old=m2(i)+1
         enddo
         do i=i1_old,nseq1
            ii=ii+1
            aseq1(ii)=seq1(i)
            aseq2(ii)='-'
            aseq3(ii)=' '
         enddo
         do i=i2_old,nseq2
            ii=ii+1
            aseq1(ii)='-'
            aseq2(ii)=seq2(i)
            aseq3(ii)=' '
         enddo
         write(*,50)
 50      format('(":" denotes the residue pairs of distance < 5.0 ',
     &        'Angstrom)')
         write(*,10)(aseq1(i),i=1,ii)
         write(*,10)(aseq3(i),i=1,ii)
         write(*,10)(aseq2(i),i=1,ii)
 10      format(10000A1)
         write(*,*)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      else if(opt%verbo_flag.eq.2) then  ! print detailed alignment

         call init_ali_map( is_ali, nseq1, nseq2 )
         do k=1,n8_al
            i=m1(k)
            j=m2(k)
            is_ali(i,j) = 1
         enddo

************  output aligned sequences in detail **************************
         write(*,*)'****** Match List ******'
!         write(*,*)'Index Ch1 Resid1  AA1 Ch2 Resid2  AA2 ',
!     &        'Distance NAC NC1 NC2'
         write(*,*)'Index Ch1 Resid1  AA1 Ch2 Resid2  AA2 ',
     &        'Distance Cos(theta)'
         tms=0
         tms1=0
         do i=1,n8_al
            dis2=sqrt((xtm1(i)-xtm2(i))**2+
     &           (ytm1(i)-ytm2(i))**2+(ztm1(i)-ztm2(i))**2)

            call get_fcol_4s(m1(i),m2(i),is_ali,fcol,col)

            i1=m1(i)
            i2=m2(i)
            iaa1=iseq1(i1)
            iaa2=iseq2(i2)

            note(1:3) = '   '
            if(seq1(i1).eq.seq2(i2).and.seq1(i1).ne.'X')then
               note(3:3) = '*'
            else if(iaa1.ge.0.and.iaa2.ge.0.and.seqmat(iaa1,iaa2).gt.0)then
               note(3:3) = ':'
            endif

            product=1
            if(iaa1.gt.0.and.iaa2.gt.0)then
               product=xtb(i1)*cbvec(1,i2,1)
               product=product+ytb(i1)*cbvec(2,i2,1)
               product=product+ztb(i1)*cbvec(3,i2,1)
            else if(iaa1.eq.0.and.iaa2.eq.0)then
               product=1        !both are glycine
            else
               product=0.272    !mean cos(theta) from random alignment
            endif
            prod=product+0.5
            if(prod.gt.1)then
               prod=1
            else if(prod.lt.0.1)then
               prod=0.1
            endif

            write(*,21) i,chain1(i1),mm1(i1),icode1(i1),seq1(i1),
     &           chain2(i2),mm2(i2),icode2(i2),seq2(i2),dis2,
     &           product,note
!            write(*,20) i,chain1(i1),mm1(i1),icode1(i1),seq1(i1),
!     &           chain2(i2),mm2(i2),icode2(i2),seq2(i2),dis2,col,
!     &           ncont1(i1),ncont2(i2),note,product
            tms=tms+1/(1+(dis2/d0)**2)
            tms1=tms1+fcol/(1+(dis2/d0)**2)
         enddo
         write(*,*)
         tms=tms/anseq
         tms1=tms1/anseq
!      tms2=tms1*op
!      write(*,'(a,f5.3,a,f5.3,a,f5.3)')'TMS = ',tms,
!     &     ', TMS (cont) = ',tms1,', TMS (order) = ',tms2
 20      format(i5,2x,2(a3,1x,i4,a1,2x,a3,2x),f6.3,i5,2i4,a4,f7.3)
 21      format(i5,2x,2(a3,1x,i4,a1,2x,a3,2x),f6.3,f11.3,a4)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      endif

      write(numstr, '(I)')nint(anseq)
      write(*,109)trim(adjustl(numstr)),d0
 109  format('Scoring parameters: normalization length = ',A,',  d0 =',F6.3)

      if(opt%order_flag.eq.0) then
         write(*,'(A)')'Alignment search mode: non-sequential'
      else
         write(*,'(A)')'Alignment search mode: sequential'
      endif

      write(*,110)best_init
 110  format('Best alignment search: initial =',I3)

      write(numstr, '(I)')tot_seqsc
      write(numstr2,'(I)')tot_seqsc_pos
      write(*,112)trim(adjustl(numstr)),trim(adjustl(numstr2))
 112  format('Sequence similarity (BLOSUM62): sum = ',A,', sum_pos = ',A)
      call system_clock ( clock_stop, clock_rate, clock_max )
      total_time = (clock_stop - clock_start)/real(clock_rate, kind=8)

      write(numstr, '(F12.5)')total_time
      write(*,'(3A,/)')'Running time: ', trim(adjustl(numstr)), ' seconds'

      RETURN
 9999 END



***********************************************************************
***********************************************************************
*     Optimal structure alignment
***********************************************************************
***********************************************************************
***********************************************************************
      SUBROUTINE super_align(invmap0,score_flag,order_flag,quick_flag,
     &     best_init,best_mode,prealn)
      use Constants
      use AlnResult

      type(alignment) :: prealn

      COMMON/BACKBONE/XA(3,maxr,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(maxr,maxr),gap_open
      common/TM/TM,TMmax
      common/contlst/cont1,cont2
      common/contnum/ncont1(maxr),ncont2(maxr)      !number of contacts

      real    coor1(3,maxr),coor2(3,maxr)
      real    gapp(100),gapp45(10)
      real    TMcut1,TMcut2  !cutoff for determing iterative runs
      real    TM_old,xa
      real    GL_i,GL_ii,GL_ir,GL_iir

      integer cont1,cont2,ncont1,ncont2
      integer nseq1,nseq2,nseq
      integer invmap0(maxr),invmap_i(maxr),invmap_ii(maxr)
      integer invmap_ir(maxr),invmap_iir(maxr)
      integer score_flag, order_flag,quick_flag,sec_flag
      integer mode,best_mode
      integer best_init
      integer n_gapp,n_gapp45,i_gapp,id,niter
      integer n_frag

      TMmax=0
      TM_old=0
      n_gapp=2
      gapp(1)=-0.6
      gapp(2)=0
      niter=20

      gap_open=0

      gapp45(1)=0
      n_gapp45=1

      TMcut1=-0.1
      TMcut2=0.2

      best_init=-1
      best_mode=0

      nseq=max(nseq1,nseq2)

      !!!coordinates of two structures
      do i=1,nseq1
         do j=1,3
            coor1(j,i)=xa(j,i,0)
         enddo
      enddo

      do i=1,nseq2
         do j=1,3
            coor2(j,i)=xa(j,i,1)
         enddo
      enddo



      !if(score_flag.eq.2) goto 10   ! for testing purpose


*000000000000000000000000000000000000000000000000000000000
*     use pre-defined alignment as the initial seeds
**********************************************************
      do 1 i=1,1
         if(prealn%length.lt.3) goto 1   ! ignore short alignment with < 3

         do j=1,nseq2
            invmap_i(j)=-1
         enddo

         do j=1,prealn%length
            m1=prealn%map1(j)
            m2=prealn%map2(j)
            invmap_i(m2)=m1
         enddo

         TM_old=TMmax
         call iter_search(TMmax,TMcut1,score_flag,quick_flag,niter,
     &        gapp,n_gapp,invmap0,invmap_i,nseq,nseq1,nseq2,0)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=0
            TM_old=TMmax
         endif
         !print '(A10,F8.5)','TMmax =',TMmax
        !call print_invmap(invmap_i,nseq2,1)

         if(order_flag.eq.0)then
            call iter_lsap(score_flag,quick_flag,invmap_i,invmap0,
     &           score,xa,nseq1,nseq2,3,mode_flag,best_mode)
            if(TMmax-TM_old.gt.1e-6)then
               best_init=0
               TM_old=TMmax
            endif
         endif
         !print '(A10,F8.5)','Initial  0:  TMmax =',TMmax

 1    enddo


*11111111111111111111111111111111111111111111111111111111
*     get initial alignment from gapless threading
*     run twice (normal and flip) to avoid asymetric score
**********************************************************
      call get_initial(score_flag,nseq1,nseq2,invmap_i,invmap_ii,
     &     coor1,coor2,cont1,cont2,ncont1,ncont2,GL_i,GL_ii)

      !!! flip order of struct 1 and struct2 to avoid asymetric score
      call get_initial(score_flag,nseq2,nseq1,invmap_ir,invmap_iir,
     &     coor2,coor1,cont2,cont1,ncont2,ncont1,GL_ir,GL_iir)
      !print '(4(A,F8.3))','GL_ii=',GL_ii,', GL_iir=',GL_iir,', GL_i=',GL_i,', GL_ir=',GL_ir

      if(GL_ii.lt.GL_iir)then
         call rev_invmap(invmap_iir,invmap_ii,nseq1,nseq2)
      endif

      if(GL_i.lt.GL_ir)then
         call rev_invmap(invmap_ir,invmap_i,nseq1,nseq2)
      endif

      !call print_invmap(invmap_ii,nseq2,1)

      TM_old=TMmax
      call iter_search(TMmax,TMcut1,score_flag,quick_flag,niter,
     &     gapp,n_gapp,invmap0,invmap_ii,nseq,nseq1,nseq2,12)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=12
         TM_old=TMmax
      endif

      call iter_search(TMmax,TMcut2,score_flag,quick_flag,niter,
     &     gapp,n_gapp,invmap0,invmap_i,nseq,nseq1,nseq2,11)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=11
         TM_old=TMmax
      endif

      !!! iterative non-sequential alignment
      if(order_flag.eq.0)then
         call iter_lsap(score_flag,quick_flag,invmap_ii,invmap0,
     &        score,xa,nseq1,nseq2,3,mode_flag,best_mode)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=12
            TM_old=TMmax
         endif

         call iter_lsap(score_flag,quick_flag,invmap_i,invmap0,
     &        score,xa,nseq1,nseq2,3,mode_flag,best_mode)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=11
            TM_old=TMmax
         endif
      endif

      !call print_invmap(invmap0,nseq2,1)


*222222222222222222222222222222222222222222222222222222222222222
*     get initial alignment from secondary structure alignment
****************************************************************
      call get_initial2(invmap_i)         !DP for secondary structure

      !call print_invmap(invmap_i,nseq2,1)
      call iter_search(TMmax,TMcut2,score_flag,quick_flag,niter,
     &     gapp,n_gapp,invmap0,invmap_i,nseq,nseq1,nseq2,2)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=2
         TM_old=TMmax
      endif

      if(order_flag.eq.0)then
         call iter_lsap(score_flag,quick_flag,invmap_i,invmap0,
     &        score,xa,nseq1,nseq2,3,mode_flag,best_mode)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=2
            TM_old=TMmax
         endif
      endif


*3333333333333333333333333333333333333333333333333333333333333
*     fragment-based alignment
**************************************************************
      n_frag=-1
      call get_initial3(invmap_i,coor1,coor2,nseq1,nseq2,n_frag,score_flag)

      call iter_search(TMmax,TMcut1,score_flag,quick_flag,niter,
     &     gapp45,n_gapp45,invmap0,invmap_i,nseq,nseq1,nseq2,32)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=32
         TM_old=TMmax
      endif

      if(order_flag.eq.0)then
         call iter_lsap(score_flag,quick_flag,invmap_i,invmap0,
     &        score,xa,nseq1,nseq2,32,mode_flag,best_mode)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=32
            TM_old=TMmax
         endif
      endif
      !call print_invmap(invmap0,nseq2,1)

      n_frag=20
      call get_initial3(invmap_i,coor1,coor2,nseq1,nseq2,n_frag,score_flag)
      !call print_invmap(invmap_i,nseq2,1)

      call iter_search(TMmax,TMcut1,score_flag,quick_flag,niter,
     &     gapp45,n_gapp45,invmap0,invmap_i,nseq,nseq1,nseq2,31)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=31
         TM_old=TMmax
      endif

      if(order_flag.eq.0)then
         call iter_lsap(score_flag,quick_flag,invmap_i,invmap0,
     &        score,xa,nseq1,nseq2,31,mode_flag,best_mode)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=31
            TM_old=TMmax
         endif
      endif

 3    continue

*444444444444444444444444444444444444444444444444444444444444
*     get initial alignment from invmap0+SS
*************************************************************
      call get_initial4(invmap_i)         !invmap0+SS

      call iter_search(TMmax,TMcut2,score_flag,quick_flag,niter,
     &     gapp,n_gapp,invmap0,invmap_i,nseq,nseq1,nseq2,4)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=4
         TM_old=TMmax
      endif

      if(order_flag.eq.0)then
         call iter_lsap(score_flag,quick_flag,invmap_i,invmap0,
     &        score,xa,nseq1,nseq2,3,mode_flag,best_mode)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=4
            TM_old=TMmax
         endif
      endif


*5555555555555555555555555555555555555555555555555555555555555555
*     get initial alignment from gapless threading with PS-score
*****************************************************************
      if( score_flag.ne.2 ) goto 6  ! only for the PS-score
      call get_initial5(score_flag,nseq1,nseq2,invmap_i,GL_i)

      call iter_search(TMmax,TMcut2,score_flag,quick_flag,niter,
     &     gapp,n_gapp,invmap0,invmap_i,nseq,nseq1,nseq2,5)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=5
         TM_old=TMmax
      endif

      if(order_flag.eq.0)then
         call iter_lsap(score_flag,quick_flag,invmap_i,invmap0,
     &        score,xa,nseq1,nseq2,3,mode_flag,best_mode)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=5
            TM_old=TMmax
         endif
      endif


      !write(*,'(a,f9.5,a,i3)')'Best sequential alignment: TMmax=',TMmax,' init=',best_init

 6    continue

*****************************************************************
*       initerative sequence order independent alignment
*****************************************************************

      if( order_flag.ne.0 ) goto 14  ! skip sequence order independent part

 10   continue


*6666666666666666666666666666666666666666666666666666666666666666
*     get initial alignment from local contact pattern
*****************************************************************
      call get_initial6(score_flag,invmap_i,invmap_ii)

      !call print_invmap(invmap_i, nseq2,1)
      call iter_lsap(score_flag,quick_flag,invmap_i,invmap0,
     &     score,xa,nseq1,nseq2,3,mode_flag,best_mode)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=61
         TM_old=TMmax
      endif

      !call print_invmap(invmap_ii,nseq2,1)
      call iter_lsap(score_flag,quick_flag,invmap_ii,invmap0,
     &     score,xa,nseq1,nseq2,3,mode_flag,best_mode)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=62
         TM_old=TMmax
      endif


 12   continue

*7777777777777777777777777777777777777777777777777777777777777777
*     get initial alignment from local contact pattern
*****************************************************************
      call get_initial7(score_flag,invmap_i,invmap_ii)
      !call print_invmap(invmap_ii,nseq2,1)

      call iter_lsap(score_flag,quick_flag,invmap_i,invmap0,
     &     score,xa,nseq1,nseq2,3,mode_flag,best_mode)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=71
         TM_old=TMmax
      endif

      call iter_lsap(score_flag,quick_flag,invmap_ii,invmap0,
     &     score,xa,nseq1,nseq2,3,mode_flag,best_mode)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=72
         TM_old=TMmax
      endif


ccc   start from the best alignment found so far
      do j=1,nseq2
         invmap_i(j)=invmap0(j)
      enddo
      !call print_invmap(invmap_i,nseq2,0)

      call iter_lsap(score_flag,quick_flag,invmap_i,invmap0,
     &     score,xa,nseq1,nseq2,30,mode_flag,best_mode)


 14   continue
c^^^^^^^^^^^^^^^ best alignment invmap0(j) found ^^^^^^^^^^^^^^^^^^
      RETURN
      END




cccccc======================================================ccccc
cccccc   Initerative dynamics programming sequential search
cccccc======================================================ccccc
      subroutine iter_search(TMmax,TMcut,score_flag,quick_flag,niter,
     &  gapp,n_gapp,invmap0,invmap_i,nseq,nseq1,nseq2,init)

      implicit none
      include  'pars.h'

      real     TM,TMmax,TM_old(10),rTMmax
      real     TMcut
      real     gapp(n_gapp),gap_open
      real     score(maxr,maxr),score0(maxr,maxr)
      real     diff,eps
      integer  niter,n_gapp,i_gapp
      integer  nseq1,nseq2,nseq
      integer  np,init
      integer  invmap(maxr),invmap0(maxr),invmap_i(maxr)
      integer  score_flag
      integer  i,j,id,ii,jj
      integer  mode,mode_flag,quick_flag

      logical  score_upd

      eps=1.E-6
      mode_flag=0
      gap_open=0

      do i=1,nseq2
         invmap(i)=invmap_i(i)
      enddo

      !call print_invmap(invmap,nseq2,1)
      call get_score(invmap,TM,score0,score_flag,mode_flag,
     &     quick_flag,mode)
      !write(*,'(a,i4,a,2F11.7)')'Initial ',init,' before DP',TM,TMmax


      if(TM - TMmax.gt.eps)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
         do j=1,nseq2
            invmap_i(j)=invmap(j)
         enddo
      endif


      !for checking convergence and periodic cases, consider a period up to 5
      np=5
      do i=1,np
         TM_old(i)=-i
      enddo

      !call print_invmap(invmap,nseq2)

      score_upd = .true.
      !Iterative dynamic programming
      if(TM-TMmax*TMcut.gt.eps.or.TMcut.lt.0)then
         DO 111 i_gapp=1,n_gapp	!different gap panalties
            GAP_OPEN=gapp(i_gapp) !gap open panalty, no gap extension penalty here
	    rTMmax=-1

            if(score_upd)then
               do ii=1,nseq1
                  do jj=1,nseq2
                     score(ii,jj)=score0(ii,jj)
                  enddo
               enddo
               score_upd=.false.
            endif

            do 222 id=1,niter

               !dynamic programming. the alignment obtained is
               !saved in invmap(j), second seq -> first seq
               call DP(score,gap_open,nseq1,nseq2,invmap)


               !calculate TM-score/IS-score, score(i,j)
               call get_score(invmap,TM,score,score_flag,
     &              mode_flag,quick_flag,mode)
               !write(*,'(I4,F6.1,I4,2F11.7)')init,gapp(i_gapp),id,TM,TMmax


               if(TM-TMmax.gt.eps)then
                  TMmax=TM
                  do j=1,nseq2
                     invmap0(j)=invmap(j)
                  enddo
                  do j=1,nseq2
                     invmap_i(j)=invmap(j)
                  enddo
               endif

	       if(TM.gt.rTMmax)then
                  rTMmax=TM
                  do ii=1,nseq1
                     do jj=1,nseq2
                        score0(ii,jj)=score(ii,jj)
                     enddo
                  enddo
                  score_upd=.true.
     	       endif

               !check convergence and periodic cases
               do i=1,np
                  if(id.gt.i.and.abs(TM-TM_old(i)).lt.eps)goto 111
               enddo

               do i=np,2,-1
                  TM_old(i)=TM_old(i-1)
               enddo
               TM_old(1)=TM

 222        continue
 111     continue
      endif
      !write(*,'(a,i4,a,F11.7)')'Initial ',init,'  after DP: TMmax=',TMmax

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccc======================================================ccccc
cccccc   Initerative non-sequential alignment search
cccccc======================================================ccccc
      subroutine iter_lsap(score_flag,quick_flag,invmap,invmap0,
     &     score,xa,nseq1,nseq2,nround,mode_flag,best_mode)
      implicit none
      include 'pars.h'

      common/TM/TM,TMmax

      real score(maxr,maxr),gap_open,z
      real xa(3,maxr,0:1)
      real TM,TMmax,TM_old(10)

      integer invmap(maxr),invmap0(maxr)
      integer nseq1,nseq2,n,naln
      integer score_flag,quick_flag,mode_flag
      integer nround,mode,best_mode
      integer id,np,i,j

      !call print_invmap(invmap,nseq2,1)

      mode_flag=1
      n=max(nseq1,nseq2)

      !call print_invmap(invmap,nseq2,1)
      call get_score(invmap,TM,score,score_flag,mode_flag,quick_flag,mode)


      !write(*,'(2(a,f8.5),a,I2)')'before iter_lsap: TM =',TM,' TMmax =',TMmax,' mode=',mode
      if(TM.gt.TMmax)then
         best_mode=mode
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif

      !for checking convergence and periodic cases
      np=10 !period up to 10
      do i=1,np
         TM_old(i)=-i
      enddo

      GAP_OPEN=0
      do 17 id=1,nround             !maximum interation is 200

         call solvLSAP(score,invmap,nseq1,nseq2,n,z) !produce alignment invmap(j)
         call remove_distant_pairs(invmap,nseq1,nseq2,xa)

         naln=0
         do i=1,nseq2
            if(invmap(i).gt.0)naln=naln+1
         enddo
         if(naln.lt.3)goto 13  ! alignment must have at least three residues

         call get_score(invmap,TM,score,score_flag,mode_flag,quick_flag,mode)
         !write(*,'(a,i2,2(a,f8.5))')'round = ',id,', TM =',TM,', TMamx =',TMmax

         if(TM-TMmax.gt.1e-6)then
            TMmax=TM
            best_mode=mode
            do j=1,nseq2
               invmap0(j)=invmap(j)
            enddo
         endif

         !check convergence and periodic cases
         do i=1,np
            if(id.gt.i.and.abs(TM-TM_old(i)).lt.1e-6)goto 13
         enddo

         do i=np,2,-1
            TM_old(i)=TM_old(i-1)
         enddo
         TM_old(1)=TM

 17   continue
 13   continue

      !write(*,'(2(a,f8.5),a,I2)')'after  iter_lsap: TM =',TM,' TMmax =',TMmax,' mode=',mode
      return
      end


**************************************************************
*     get initial alignment from gapless threading
**************************************************************
      subroutine get_initial(score_flag,nseq1,nseq2,invmap_i,invmap_ii,
     &     coor1,coor2,cont1,cont2,ncont1,ncont2,GL_max,GL_maxA)
      implicit none
      include 'pars.h'

      common/d0/d0,anseq
      common/d0min/d0_min
      common/dpc/score(maxr,maxr),gap_open

      real    coor1(3,maxr),coor2(3,maxr)
      real    score,gap_open
      integer nseq1,nseq2,nseq   !length of full structures
      integer invmap(maxr),invmap_i(maxr),invmap_ii(maxr)
      integer cont1(maxc,maxr),cont2(maxc,maxr)
      integer ncont1(maxr),ncont2(maxr)

      real    d0,anseq,d0_min
      real    d01,d02,xx,yy,zz,dd
      real    aL,GL_cf,GL,GL_max,GL_maxA


      double precision r_1(3,maxr),r_2(3,maxr),r_3(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms,drms
      data w /maxr*1.0/
      integer ier

      real  xtmm1(maxr),ytmm1(maxr),ztmm1(maxr)
      real  xtmm2(maxr),ytmm2(maxr),ztmm2(maxr)

      integer n_jump,n1,n2,idel,ishift
      integer i,j,L,ii,i1,j1
      integer score_flag


      d01=d0+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01*d01

      n_jump=1
      aL=min(nseq1,nseq2)
      idel=nint(aL/2.0)             !minimum size of considered fragment
      if(idel.le.5)idel=5
      n1=-nseq2+idel
      n2=nseq1-idel
      GL_max=0
      GL_maxA=0
      GL_cf=idel/2.
      gap_open=0.

      do j=1,nseq2
         invmap(j) = -1
         invmap_i(j) = -1
         invmap_ii(j)= -1
      enddo

      do ishift=n1,n2,n_jump
         L=0
         do j=1,nseq2
            i=j+ishift
            if(i.ge.1.and.i.le.nseq1)then
               L=L+1
               invmap(j)=i
               xtmm1(L)=coor1(1,i)
               ytmm1(L)=coor1(2,i)
               ztmm1(L)=coor1(3,i)
               xtmm2(L)=coor2(1,j)
               ytmm2(L)=coor2(2,j)
               ztmm2(L)=coor2(3,j)
            else
               invmap(j)=-1
            endif
         enddo

         if(L.ge.idel)then
            call get_GL(d0,GL,nseq1,nseq2,invmap,coor1,coor2,
     &     cont1,cont2,ncont1,ncont2,0)
            if(GL.gt.GL_max*0.95)then
               !print '(A,F8.3,A,F8.3,A,I4)','GL=',GL,', GL_max=',GL_max, ', nseq2=',nseq2
               if(GL.gt.GL_max)then
                  GL_max=GL
                  do i=1,nseq2
                     invmap_i(i)=invmap(i)
                  enddo
                  !call print_invmap(invmap_i,nseq2,1)
               endif

               do ii=1,L
                  r_1(1,ii)=xtmm1(ii)
                  r_1(2,ii)=ytmm1(ii)
                  r_1(3,ii)=ztmm1(ii)
                  r_2(1,ii)=xtmm2(ii)
                  r_2(2,ii)=ytmm2(ii)
                  r_2(3,ii)=ztmm2(ii)
               enddo

               call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2

               do i1=1,nseq1
                  xx=t(1)+u(1,1)*coor1(1,i1)+u(1,2)*coor1(2,i1)+u(1,3)*
     &                 coor1(3,i1)
                  yy=t(2)+u(2,1)*coor1(1,i1)+u(2,2)*coor1(2,i1)+u(2,3)*
     &                 coor1(3,i1)
                  zz=t(3)+u(3,1)*coor1(1,i1)+u(3,2)*coor1(2,i1)+u(3,3)*
     &                 coor1(3,i1)
                  do j1=1,nseq2
                     dd=(xx-coor2(1,j1))**2+(yy-coor2(2,j1))**2+(zz-
     &                    coor2(3,j1))**2
                     score(i1,j1)=1/(1+dd/d02)
                  enddo
               enddo

*********extract alignement with score(i,j) *****************

               call DP(score,gap_open,NSEQ1,NSEQ2,invmap)
               call get_GL(d0,GL,nseq1,nseq2,invmap,coor1,coor2,
     &              cont1,cont2,ncont1,ncont2,0)

               if(GL.gt.GL_maxA)then
                  !print '(A,F8.3,A,F8.3,A,I4)','GL=',GL,', GL_maxA=',GL_maxA, ', nseq2=',nseq2
                  GL_maxA=GL
                  do i1=1,nseq2
                     invmap_ii(i1)=invmap(i1)
                  enddo
                  !call print_invmap(invmap_ii,nseq2,1)
               endif

            endif
         endif
      enddo

      !print *,'At the end of get_intial: '
      !call print_invmap(invmap_i,nseq2,1)

      return
      end

**************************************************************
*     get initial alignment invmap0(i) from secondary structure
**************************************************************
      subroutine get_initial2(invmap)
      include 'pars.h'
      COMMON/BACKBONE/XA(3,maxr,0:1)
      common/length/nseq1,nseq2
      common/alignrst/invmap0(maxr)
      common/TM/TM,TMmax
      common/sec/isec(maxr),jsec(maxr)
      common/dpc/score(maxr,maxr),gap_open

      integer invmap(maxr)


********** score matrix **************************
      do i=1,nseq2
         invmap(i)=invmap0(i)
      enddo
      call get_score1(invmap)           !get score(i,j) using RMSD martix

      do i=1,nseq1
         do j=1,nseq2
            if(isec(i).eq.jsec(j))then
               score(i,j)=1+0.01*score(i,j)
            else
               score(i,j)=0.01*score(i,j)
            endif
         enddo
      enddo

********** find initial alignment: invmap(j) ************
      gap_open=-1.0             !should be -1
      call DP(score,gap_open,NSEQ1,NSEQ2,invmap)      !produce alignment invmap(j)

*^^^^^^^^^^^^ initial alignment done ^^^^^^^^^^^^^^^^^^^^^^
      return
      end




**************************************************************
*    Fragment-based superposition.                           *
**************************************************************
      subroutine get_initial3(invmap_i,coor1,coor2,nseq1,nseq2,n_frag,score_flag)
      implicit none
      include 'pars.h'

      common/dpc/score(maxr,maxr),gap_open
      common/d0/d0,anseq
      common/d0min/d0_min
      common/sec/isec(maxr),jsec(maxr)
      common/contlst/cont1,cont2
      common/contnum/ncont1(maxr),ncont2(maxr)      !number of contacts
      integer cont1(maxc,maxr),cont2(maxc,maxr)
      integer ncont1,ncont2

      real    xa,score,gap_open
      real    d0,anseq,d0_min
      real    d01,d02
      real    coor1(3,maxr),coor2(3,maxr)

      double precision r_1(3,maxr),r_2(3,maxr),r_3(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms,drms
      data w /maxr*1.0/
      integer ier

      integer m1,m2,n_frag,n_frag2
      integer invmap(maxr),invmap_i(maxr)
      integer nseq1,nseq2
      integer n_jump1,n_jump2,n_jump
      integer isec,jsec
      integer secmatch
      integer i,j,ii,jj,k,kk,i1,j1
      integer minlen,maxlen
      integer score_flag

      real    GL,GL_cf,GLmaxA,GLmaxA_cf
      real    xx,yy,zz,dd


***** setting parameters ************************************
      gap_open=0
      maxlen=max(nseq1,nseq2)
      minlen=min(nseq1,nseq2)

      if(n_frag.lt.5) then
         n_frag=minlen/20       !length of the sliding window
         if(n_frag.lt.5)n_frag=5
         n_frag2=n_frag/2
         if(n_frag2.lt.2)n_frag2=2

         n_jump1=nseq1/20 - 2
         n_jump2=nseq2/20 - 2

         if(n_jump1.lt.1)n_jump1=1
         if(n_jump2.lt.1)n_jump2=1
         if(n_jump1.gt.n_frag)n_jump1=n_frag
         if(n_jump2.gt.n_frag)n_jump2=n_frag

         GL_cf=n_frag*0.75
      else
         if(n_frag.gt.minlen)n_frag=minlen
         n_frag2=0
         n_jump1=n_frag
         n_jump2=n_frag
         GL_cf=n_frag*0.1
      endif


      !n_jump=min(n_jump1,n_jump2)
      !print *,'n_frag=',n_frag,' n_jump1=',n_jump1,' n_jump2=',n_jump2

      m1=nseq1-n_frag+1
      m2=nseq2-n_frag+1

      d01=d0+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01*d01

      GLmaxA=0
      GL=0

      GLmaxA_cf=min(nseq1,nseq2)*0.9

      do j1=1,nseq2
         invmap_i(j1)=-1
      enddo

      do 10 ii=1,m1,n_jump1
         do 20 jj=1,m2,n_jump2
            secmatch=0
            do j=1,nseq2
               invmap(j)=-1
            enddo

            do 100 kk=1,n_frag
               i=ii+kk-1
               j=jj+kk-1
               if(i.gt.nseq1.or.j.gt.nseq2) goto 100
               invmap(j)=i
               if(isec(i).eq.jsec(j))then
                  secmatch=secmatch+1
               endif
 100        continue

            !print *,'secmatch=',secmatch,', n_frag2=',n_frag2
            if(secmatch.le.n_frag2)then
               goto 20
            endif

            !!! first identify similar fragements
            call get_GL(d0,GL,nseq1,nseq2,invmap,coor1,coor2,
     &     cont1,cont2,ncont1,ncont2,0)
            if(GL.le.GL_cf)goto 20

            !!! use good fragment alignment to align whole structures
            k=0
            do 200 j=1,nseq2
               i=invmap(j)
               if(i.le.0)goto 200
               k=k+1
               r_1(1,k)=coor1(1,i)
               r_1(2,k)=coor1(2,i)
               r_1(3,k)=coor1(3,i)
               r_2(1,k)=coor2(1,j)
               r_2(2,k)=coor2(2,j)
               r_2(3,k)=coor2(3,j)
 200        continue

*********superpose the two structures and rotate it *****************
            call u3b(w,r_1,r_2,k,1,rms,u,t,ier) !u rotate r_1 to r_2

            !call print_TransMatrix(t,u)

            do i1=1,nseq1
               xx=t(1)+u(1,1)*coor1(1,i1)+u(1,2)*coor1(2,i1)+u(1,3)
     &              *coor1(3,i1)
               yy=t(2)+u(2,1)*coor1(1,i1)+u(2,2)*coor1(2,i1)+u(2,3)
     &              *coor1(3,i1)
               zz=t(3)+u(3,1)*coor1(1,i1)+u(3,2)*coor1(2,i1)+u(3,3)
     &              *coor1(3,i1)
               do j1=1,nseq2
                  dd=(xx-coor2(1,j1))**2+(yy-coor2(2,j1))**2+
     &                 (zz-coor2(3,j1))**2
                  score(i1,j1)=1/(1+dd/d02) ! changing
               enddo
            enddo


*********extract alignement with score(i,j) *****************
            call DP(score,gap_open,NSEQ1,NSEQ2,invmap)
            call get_GL(d0,GL,nseq1,nseq2,invmap,coor1,coor2,
     &     cont1,cont2,ncont1,ncont2,0)
            !print *,'GL=',GL,' GLmaxA=',GLmaxA

            if(GL.gt.GLmaxA)then
               GLmaxA=GL
               do j1=1,nseq2
                  invmap_i(j1)=invmap(j1)
               enddo
               if(GLmaxA.gt.GLmaxA_cf)return !reduce scans for highly similar structures
            endif

 20      continue
 10   continue

      return
      end


**************************************************************
*     get initial alignment invmap0(i) from secondary structure
*     and previous alignments
**************************************************************
      subroutine get_initial4(invmap)
      include 'pars.h'
      COMMON/BACKBONE/XA(3,maxr,0:1)
      common/length/nseq1,nseq2
      common/alignrst/invmap0(maxr)
      common/TM/TM,TMmax
      common/sec/isec(maxr),jsec(maxr)
      common/dpc/score(maxr,maxr),gap_open

      integer invmap(maxr)

********** score matrix **************************
      do i=1,nseq2
         invmap(i)=invmap0(i)
      enddo
      call get_score1(invmap)           !get score(i,j) using RMSD martix
      do i=1,nseq1
         do j=1,nseq2
            if(isec(i).eq.jsec(j))then
               score(i,j)=0.5+score(i,j)
            else
               score(i,j)=score(i,j)
            endif
         enddo
      enddo

********** find initial alignment: invmap(j) ************
      gap_open=-1.0             !should be -1
      call DP(score,gap_open,NSEQ1,NSEQ2,invmap)      !produce alignment invmap(j)


*^^^^^^^^^^^^ initial alignment done ^^^^^^^^^^^^^^^^^^^^^^
      return
      end

**************************************************************
*     determine secondary structure for each residue
**************************************************************
      subroutine calcSecStruct( pdb )
      use PDBType

      type(pdbstruct) :: pdb
      real coor(3,maxr)
      integer isec(maxr)

      nseq1 = pdb%num_residue

      do i=1,nseq1
         do j=1,3
            coor(j,i) = pdb%coor_ca(j,i)
         enddo
      enddo

********** assign secondary structures ***************
c     1->coil, 2->helix, 3->turn, 4->strand
      do i=1,nseq1
         isec(i)=1
         j1=i-2
         j2=i-1
         j3=i
         j4=i+1
         j5=i+2
         if(j1.ge.1.and.j5.le.nseq1)then
            dis13=diszy(coor,j1,j3)
            dis14=diszy(coor,j1,j4)
            dis15=diszy(coor,j1,j5)
            dis24=diszy(coor,j2,j4)
            dis25=diszy(coor,j2,j5)
            dis35=diszy(coor,j3,j5)
            isec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
         endif
      enddo


***   smooth single -------------->
***   --x-- => -----
      do i=1,nseq1
         if(isec(i).eq.2.or.isec(i).eq.4)then
            j=isec(i)
            if(isec(i-2).ne.j)then
               if(isec(i-1).ne.j)then
                  if(isec(i+1).ne.j)then
                     if(isec(i+2).ne.j)then  ! Mu: i+2 instead of i+1
                        isec(i)=1
                     endif
                  endif
               endif
            endif
         endif
      enddo

***   smooth double -------------->
***   --xx-- => ------
      do i=1,nseq1
         if(isec(i).ne.2)then
         if(isec(i+1).ne.2)then
         if(isec(i+2).eq.2)then
         if(isec(i+3).eq.2)then
         if(isec(i+4).ne.2)then
         if(isec(i+5).ne.2)then
            isec(i+2)=1
            isec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif

         if(isec(i).ne.4)then
         if(isec(i+1).ne.4)then
         if(isec(i+2).eq.4)then
         if(isec(i+3).eq.4)then
         if(isec(i+4).ne.4)then
         if(isec(i+5).ne.4)then
            isec(i+2)=1
            isec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif
      enddo


***   connect -------------->
***   x-x => xxx
      do i=1,nseq1
         if(isec(i).eq.2)then
         if(isec(i+1).ne.2)then
         if(isec(i+2).eq.2)then
            isec(i+1)=2
         endif
         endif
         endif

         if(isec(i).eq.4)then
         if(isec(i+1).ne.4)then
         if(isec(i+2).eq.4)then
            isec(i+1)=4
         endif
         endif
         endif
      enddo

      do i=1,nseq1
         pdb%sec_struct(i)=isec(i)
      enddo

      return
      end

*************************************************************
*     assign secondary structure:
*************************************************************
      function diszy(coor,i1,i2)
      include 'pars.h'

      real coor(3,maxr)

      diszy=sqrt((coor(1,i1)-coor(1,i2))**2
     &     +(coor(2,i1)-coor(2,i2))**2
     &     +(coor(3,i1)-coor(3,i2))**2)
      return
      end

*************************************************************
*     assign secondary structure:
*************************************************************
      function make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
      make_sec=1
      delta=2.1
      if(abs(dis15-6.37).lt.delta)then
         if(abs(dis14-5.18).lt.delta)then
            if(abs(dis25-5.18).lt.delta)then
               if(abs(dis13-5.45).lt.delta)then
                  if(abs(dis24-5.45).lt.delta)then
                     if(abs(dis35-5.45).lt.delta)then
                        make_sec=2 !helix
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      delta=1.42
      if(abs(dis15-13).lt.delta)then
         if(abs(dis14-10.4).lt.delta)then
            if(abs(dis25-10.4).lt.delta)then
               if(abs(dis13-6.1).lt.delta)then
                  if(abs(dis24-6.1).lt.delta)then
                     if(abs(dis35-6.1).lt.delta)then
                        make_sec=4 !strand
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      if(dis15.lt.8)then
         make_sec=3
      endif

      return
      end




**************************************************************
*     get initial alignment from gapless threading
**************************************************************
      subroutine get_initial5(score_flag,nseq1,nseq2,invmap_i,
     &     GL_max)
      implicit none
      include 'pars.h'
      common/d0/d0,anseq

      integer nseq1,nseq2,nseq   !length of full structures
      integer invmap(maxr),invmap_i(maxr)

      real    aL,GL_cf,GL,GL_max
      real    d0,anseq
      real    score(maxr,maxr)

      integer n_jump,n1,n2,idel,ishift
      integer i,j,L
      integer score_flag


      n_jump=1
      aL=min(nseq1,nseq2)
      idel=nint(aL/2.0)             !minimum size of considered fragment
      if(idel.le.5)idel=5
      n1=-nseq2+idel
      n2=nseq1-idel
      GL_max=0

      do j=1,nseq2
         invmap(j) = -1
         invmap_i(j) = -1
      enddo

      do ishift=n1,n2,n_jump
         L=0
         do j=1,nseq2
            i=j+ishift
            if(i.ge.1.and.i.le.nseq1)then
               L=L+1
               invmap(j)=i
            else
               invmap(j)=-1
            endif
         enddo

         if(L.ge.idel)then
            call get_GLD(GL,3.5,score_flag,invmap)

            if(GL.gt.GL_max)then
               !print '(A,F8.3,A,F8.3,A,I4)','GL=',GL,', GL_max=',GL_max, ', nseq2=',nseq2
               if(GL.gt.GL_max)then
                  GL_max=GL
                  do i=1,nseq2
                     invmap_i(i)=invmap(i)
                  enddo
                  !call print_invmap(invmap_i,nseq2,1)
               endif
            endif
         endif
      enddo

      return
      end




**************************************************************
*     local contact alignment
**************************************************************
      subroutine get_initial6(score_flag,invmap,invmap_i)
      implicit none
      include 'pars.h'

      COMMON/BACKBONE/XA(3,maxr,0:1)
      common/length/nseq1,nseq2
      common/d0/d0,anseq
      common/cbvec/cbvec(3,maxr,0:1)  ! normalized vector from ca to cb
      common/dpc/score(maxr,maxr),gap_open
      common/contlst/cont1,cont2
      common/contnum/ncont1(maxr),ncont2(maxr)      !number of contacts

      real    xa,d0,anseq
      real    score,score1(maxr,maxr)
      real    d0_search,gap_open
      real    cbvec

      integer nseq1,nseq2,invmap(maxr),invmap_i(maxr)
      integer ncont1,ncont2
      integer isec,jsec
      integer cont1(maxc,maxr),cont2(maxc,maxr)     !contact list

      double precision r_1(3,maxr),r_2(3,maxr),r_3(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms,drms
      integer ier
      data w /maxr*1.0/

      real xtmm1(maxr),ytmm1(maxr),ztmm1(maxr)
      real xtmm2(maxr),ytmm2(maxr),ztmm2(maxr)
      real GL,GL_max,xx,yy,zz,dd
      real xxb,yyb,zzb
      real z,z_min
      real GLmax(maxr),index(maxr),product

      integer top1(maxr),top2(maxr), top_cf
      integer score_flag
      integer nGLmax
      integer i,j,k,i1,j1,len,len1,len2,nmax,lmax
      integer ii,jj,maxlen

      integer iaa1,iaa2,seqsc
      integer seqmat(0:19,0:19),blmat(0:19,0:19)   ! sequence substitution matrix, e.g., BLOSUM62
      integer iseq1(maxr),iseq2(maxr)
      real    sdel   ! penalty scaling factor for unmatched amino acids
      common/seq/iseq1,iseq2,seqmat,sdel
      common/blmat/blmat


      GL_max=0
      nGLmax=0
      z_min=1000000
      nmax=max(nseq1,nseq2)

      maxlen = 5

      d0_search = 3.5
      top_cf=100      !maximum cases to consider

      do i=1,maxr
         GLmax(i)=-1
      enddo

      do 10 i=1,nseq1
         len1=ncont1(i)
         do 20 j=1,nseq2
            iaa1=iseq1(i)
            iaa2=iseq2(j)
            score1(i,j)=0.05*blmat(iaa1,iaa2)

            len2=ncont2(j)
            len=min(len1,len2)
            if(len.lt.maxlen)goto 20

            do k=1,nseq2
               invmap(k)=-1
            enddo

            do k=1,maxlen
               i1=cont1(k,i)
               j1=cont2(k,j)
               invmap(j1)=i1
            enddo

            call get_GLD(GL,3.5,score_flag,invmap)
            GL=GL/maxlen
            score1(i,j)=score1(i,j)+0.5*GL

            !print '(3(A,F11.3))','GL = ', GL, ' GLmax=',GL_max,', score1(i,j)=',score1(i,j)
            if(nGLmax.lt.maxr)then
               if(GL.gt.GL_max)GL_max=GL
               nGLmax=nGLmax+1
               GLmax(nGLmax)=GL     !record good cases
               index(nGLmax)=nGLmax
               top1(nGLmax)=i
               top2(nGLmax)=j
            endif

 20   continue
 10   continue

      !print *,'nGLmax=',nGLmax
      call ssort(GLmax,index, nGLmax, -2)
      do 30 ii=1,top_cf
         jj=index(ii)
         if(jj.le.0) goto 30

         i=top1(jj)
         j=top2(jj)

         len1=ncont1(i)
         len2=ncont2(j)
         len=min(len1,len2)
         len=min(len,maxlen)

         do k=1,len
            i1=cont1(k,i)
            j1=cont2(k,j)
            r_1(1,k)=xa(1,i1,0)
            r_1(2,k)=xa(2,i1,0)
            r_1(3,k)=xa(3,i1,0)
            r_2(1,k)=xa(1,j1,1)
            r_2(2,k)=xa(2,j1,1)
            r_2(3,k)=xa(3,j1,1)
         enddo

         call u3b(w,r_1,r_2,len,1,rms,u,t,ier) !u rotate r_1 to r_2

         do i1=1,nseq1
            xx=t(1)+u(1,1)*xa(1,i1,0)+u(1,2)*xa(2,i1,0)+u(1,3)*xa(3,i1,0)
            yy=t(2)+u(2,1)*xa(1,i1,0)+u(2,2)*xa(2,i1,0)+u(2,3)*xa(3,i1,0)
            zz=t(3)+u(3,1)*xa(1,i1,0)+u(3,2)*xa(2,i1,0)+u(3,3)*xa(3,i1,0)

            xxb=u(1,1)*cbvec(1,i1,0)+u(1,2)*cbvec(2,i1,0)+u(1,3)*cbvec(3,i1,0)
            yyb=u(2,1)*cbvec(1,i1,0)+u(2,2)*cbvec(2,i1,0)+u(2,3)*cbvec(3,i1,0)
            zzb=u(3,1)*cbvec(1,i1,0)+u(3,2)*cbvec(2,i1,0)+u(3,3)*cbvec(3,i1,0)

            do j1=1,nseq2
               dd=(xx-xa(1,j1,1))**2+(yy-xa(2,j1,1))**2+(zz-xa(3,j1,1))**2
               if(score_flag.eq.2) then
                  iaa1=iseq1(i1)
                  iaa2=iseq2(j1)
                  if(iaa1.ge.0.and.iaa2.ge.0)then
                     seqsc=seqmat(iaa1,iaa2)
                     if(seqsc.gt.0)then
                        score(i1,j1)=seqsc/(1+dd/d0_search**2)
                     else
                        score(i1,j1)=sdel/(1+dd/d0_search**2)
                     endif
                     if(iaa1.gt.0.and.iaa2.gt.0)then
                        product=xxb*cbvec(1,j1,1)
                        product=product+yyb*cbvec(2,j1,1)
                        product=product+zzb*cbvec(3,j1,1)
                     else if(iaa1.eq.0.and.iaa2.eq.0)then
                        product=1 !both are glycine
                     else
                        product=0.272   ! mean cos(theta) from random alignment
                     endif
                     product=product+0.5
                     if(product.gt.1)then
                        product=1
                     else if(product.lt.0.1)then
                        product=0.1
                     endif
                     score(i1,j1)=score(i1,j1)*product
                  else
                     score(i1,j1)=sdel/(1+dd/d0_search**2)
                  endif
               else
                  score(i1,j1)=1/(1+dd/d0_search**2)
               endif

            enddo
         enddo

         call solvLSAP(score,invmap,nseq1,nseq2,nmax,z)
         call remove_distant_pairs(invmap,nseq1,nseq2,xa)

         !write(*,'(3(A,I5),3(A,F12.4))')'i=',i,' j=',j,' len=',len,' GL=',GLmax(ii),' z=',z,' z_min=',z_min
         if(z.lt.z_min)then
            !write(*,'(A,I5,2(A,F11.2))')'lmax=',lmax,' z=',z,' z_min=',z_min
            z_min=z
            do j1=1,nseq2
               invmap_i(j1)=invmap(j1)
            enddo
            !call print_invmap(invmap_i,nseq2,1)
         endif
 30   continue

      do j1=1,nseq2
         invmap(j1)=invmap_i(j1)
      enddo
      !call print_invmap(invmap,nseq2,0)

      call get_score1(invmap)

      do 40 i=1,nseq1
         do 50 j=1,nseq2
            score1(i,j)=score1(i,j)+0.001*score(i,j)   ! make optimal alignment unique
 50      enddo
 40   enddo

      call solvLSAP(score1,invmap_i,nseq1,nseq2,nmax,z)
      !write(*,'(A,F11.5)')' z=',z
      !call print_invmap(invmap_i,nseq2,0)
      call remove_distant_pairs(invmap_i,nseq1,nseq2,xa)


      return
      end






**************************************************************
*     local contact alignment
**************************************************************
      subroutine get_initial7(score_flag,invmap,invmap_i)
      implicit none
      include 'pars.h'

      COMMON/BACKBONE/XA(3,maxr,0:1)
      common/length/nseq1,nseq2
      common/d0/d0,anseq
      common/cbvec/cbvec(3,maxr,0:1)  ! normalized vector from ca to cb
      common/dpc/score(maxr,maxr),gap_open
      common/contlst/cont1,cont2
      common/contnum/ncont1(maxr),ncont2(maxr)      !number of contacts
      common/contdis/dcont1,dcont2                  !contact distance

      real    dcont1(maxr,maxr),dcont2(maxr,maxr)
      real    xa,d0,anseq
      real    score,score1(maxr,maxr)
      real    d0_search,gap_open
      real    cbvec

      integer nseq1,nseq2,invmap(maxr),invmap_i(maxr)
      integer ncont1,ncont2
      integer cont1(maxc,maxr),cont2(maxc,maxr)     !contact list
      integer inv(100,maxr)

      double precision r_1(3,maxr),r_2(3,maxr),r_3(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms,drms
      integer ier
      data w /maxr*1.0/

      real xtmm1(maxr),ytmm1(maxr),ztmm1(maxr)
      real xtmm2(maxr),ytmm2(maxr),ztmm2(maxr)
      real GL,GL_max,xx,yy,zz,dd,dd1,dd2,prev_GL
      real xxb,yyb,zzb
      real z,z_min,prev_z
      real GLmax(maxr),index(maxr),product,GLmin

      integer top1(maxr),top2(maxr), top_cf
      integer score_flag
      integer nGLmax,GLmin_ind
      integer i,j,k,k1,k2,i1,j1,len,len1,len2,nmax,lmax
      integer ii,jj,maxlen,maxlen1,maxlen2

      integer flag
      integer iaa1,iaa2,seqsc
      integer seqmat(0:19,0:19),blmat(0:19,0:19)   ! sequence substitution matrix, e.g., BLOSUM62
      integer iseq1(maxr),iseq2(maxr)
      real    sdel   ! penalty scaling factor for unmatched amino acids
      common/seq/iseq1,iseq2,seqmat,sdel
      common/blmat/blmat


      GL_max=0
      nGLmax=0
      z_min=1000000
      nmax=max(nseq1,nseq2)

      maxlen = 7

      d0_search = 3.5
      top_cf=100      !maximum cases to consider

      do i=1,maxr
         GLmax(i)=-1
      enddo

      prev_GL=-1
      do 10 i=1,nseq1
         len1=ncont1(i)
         maxlen1=min(len1,maxlen)
         do 20 j=1,nseq2
            iaa1=iseq1(i)
            iaa2=iseq2(j)
            score1(i,j)=0.05*blmat(iaa1,iaa2)

            len2=ncont2(j)
            maxlen2=min(len2,maxlen)
            len=min(maxlen1,maxlen2)
            if(len.lt.4)goto 20

            do k=1,nseq2
               invmap(k)=-1
            enddo

            do k1=1,maxlen1
               i1=cont1(k1,i)
               dd1=dcont1(i1,i)
               do k2=1,maxlen2
                  j1=cont2(k2,j)
                  dd2=dcont2(j1,j)
                  score(k1,k2)=0.5/(1+((dd1-dd2)/d0_search)**2)
     &                 +0.05*blmat(iseq1(i1),iseq2(j1))
                  !if(i.eq.2.and.j.eq.2)print '(2(A,I4),A,F8.3)','i1=',i1,', j1=', j1,', score(k1,k2)=',score(k1,k2)
               enddo
            enddo

            len=max(maxlen1,maxlen2)
            call solvLSAP(score,invmap,maxlen1,maxlen2,len,z)

            ! convert back to the original sequence indexes
            do k2=1,nseq2
               invmap_i(k2)=-1
            enddo
            do k2=1,nseq2
               k1=invmap(k2)
               if(k1.gt.0)then
                  i1=cont1(k1,i)
                  j1=cont2(k2,j)
                  invmap_i(j1)=i1
               endif
            enddo
            call remove_distant_pairs(invmap_i,nseq1,nseq2,xa)

            GL=2-z/maxlen
            score1(i,j)=score1(i,j)+0.5*GL
            !print '(2(A,I4),A,F8.4)','i =',i,', j =',j,', GL =', GL

            if(nGLmax.lt.top_cf)then
               if(GL.gt.GL_max)GL_max=GL
               nGLmax=nGLmax+1
               GLmax(nGLmax)=GL     !record good cases
               index(nGLmax)=nGLmax
               top1(nGLmax)=i
               top2(nGLmax)=j
               do k=1,nseq2
                  inv(nGLmax,k)=invmap_i(k)
               enddo
            else !if more than top_cf, find the one with smallest GL, replace it with the current one if the current GL > smallest GL
               GLmin=1000000
               GLmin_ind=-1
               do k1=1,top_cf
                  if(GLmax(k1).lt.GLmin)then
                     GLmin=GLmax(k1)
                     GLmin_ind=k1
                  endif
               enddo
               if(GL.gt.GLmin.and.GLmin_ind.gt.0)then
                  GLmax(GLmin_ind)=GL
                  top1(GLmin_ind)=i
                  top2(GLmin_ind)=j
                  do k=1,nseq2
                     inv(GLmin_ind,k)=invmap_i(k)
                  enddo
               endif
            endif

 20   continue
 10   continue

      !print *,'nGLmax = ',nGLmax
      call ssort(GLmax,index, nGLmax, -2)

      do 30 ii=1,top_cf
         jj=index(ii)
         GL=GLmax(ii)
         if(jj.le.0) goto 30

         len=0
         do j=1,nseq2
            i=inv(jj,j)
            if(i.gt.0)then
               len=len+1
               r_1(1,len)=xa(1,i,0)
               r_1(2,len)=xa(2,i,0)
               r_1(3,len)=xa(3,i,0)
               r_2(1,len)=xa(1,j,1)
               r_2(2,len)=xa(2,j,1)
               r_2(3,len)=xa(3,j,1)
            endif
         enddo

         call u3b(w,r_1,r_2,len,1,rms,u,t,ier) !u rotate r_1 to r_2

         do i1=1,nseq1
            xx=t(1)+u(1,1)*xa(1,i1,0)+u(1,2)*xa(2,i1,0)+u(1,3)*xa(3,i1,0)
            yy=t(2)+u(2,1)*xa(1,i1,0)+u(2,2)*xa(2,i1,0)+u(2,3)*xa(3,i1,0)
            zz=t(3)+u(3,1)*xa(1,i1,0)+u(3,2)*xa(2,i1,0)+u(3,3)*xa(3,i1,0)

            xxb=u(1,1)*cbvec(1,i1,0)+u(1,2)*cbvec(2,i1,0)+u(1,3)*cbvec(3,i1,0)
            yyb=u(2,1)*cbvec(1,i1,0)+u(2,2)*cbvec(2,i1,0)+u(2,3)*cbvec(3,i1,0)
            zzb=u(3,1)*cbvec(1,i1,0)+u(3,2)*cbvec(2,i1,0)+u(3,3)*cbvec(3,i1,0)

            do j1=1,nseq2
               dd=(xx-xa(1,j1,1))**2+(yy-xa(2,j1,1))**2+(zz-xa(3,j1,1))**2
               if(score_flag.eq.2) then
                  iaa1=iseq1(i1)
                  iaa2=iseq2(j1)
                  if(iaa1.ge.0.and.iaa2.ge.0)then
                     seqsc=seqmat(iaa1,iaa2)
                     if(seqsc.gt.0)then
                        score(i1,j1)=seqsc/(1+dd/d0_search**2)
                     else
                        score(i1,j1)=sdel/(1+dd/d0_search**2)
                     endif
                     if(iaa1.gt.0.and.iaa2.gt.0)then
                        product=xxb*cbvec(1,j1,1)
                        product=product+yyb*cbvec(2,j1,1)
                        product=product+zzb*cbvec(3,j1,1)
                     else if(iaa1.eq.0.and.iaa2.eq.0)then
                        product=1 !both are glycine
                     else
                        product=0.272   ! mean cos(theta) from random alignment
                     endif
                     product=product+0.5
                     if(product.gt.1)then
                        product=1
                     else if(product.lt.0.1)then
                        product=0.1
                     endif
                     score(i1,j1)=score(i1,j1)*product
                  else
                     score(i1,j1)=sdel/(1+dd/d0_search**2)
                  endif
               else
                  score(i1,j1)=1/(1+dd/d0_search**2)
               endif

            enddo
         enddo

         call solvLSAP(score,invmap,nseq1,nseq2,nmax,z)
         call remove_distant_pairs(invmap,nseq1,nseq2,xa)

         !write(*,'(3(A,I5),3(A,F12.4))')'i=',top1(jj),' j=',top2(jj),' len=',len,' GL=',GLmax(ii),' z=',z,' z_min=',z_min
         if(z.lt.z_min)then
            !write(*,'(A,I5,2(A,F11.2))')'round =',ii,' z=',z,' z_min=',z_min
            z_min=z
            do j1=1,nseq2
               invmap_i(j1)=invmap(j1)
            enddo
            !call print_invmap(invmap_i,nseq2,1)
         endif
 30   continue

      do j1=1,nseq2
         invmap(j1)=invmap_i(j1)
      enddo
      !call print_invmap(invmap,nseq2,0)

      call get_score1(invmap)

      do 40 i=1,nseq1
         do 50 j=1,nseq2
            score1(i,j)=score1(i,j)+0.001*score(i,j)   ! make optimal alignment unique
 50      enddo
 40   enddo

      call solvLSAP(score1,invmap_i,nseq1,nseq2,nmax,z)
      !write(*,'(A,F11.5)')' z=',z
      !call print_invmap(invmap_i,nseq2,0)
      call remove_distant_pairs(invmap_i,nseq1,nseq2,xa)


      return
      end






**************************************************************
*     transform coordinates of structure 1 according to invmap
**************************************************************
      subroutine transcoor(coor,invmap,nseq1,nseq2)
      implicit none
      include 'pars.h'
      COMMON/BACKBONE/XA(3,maxr,0:1)

      real    coor(3,maxr,0:1),xa
      real    xx(maxr),yy(maxr),zz(maxr)
      integer nseq1,nseq2,invmap(maxr),nal

      double precision r_1(3,maxr),r_2(3,maxr),r_3(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms,drms
      integer ier
      data w /maxr*1.0/

      integer i,j,k,i1

      nal=0
      do 10 j=1,nseq2
         i=invmap(j)
         if(i.le.0) goto 10
         nal=nal+1
         r_1(1,nal)=xa(1,i,0)
         r_1(2,nal)=xa(2,i,0)
         r_1(3,nal)=xa(3,i,0)
         r_2(1,nal)=xa(1,j,1)
         r_2(2,nal)=xa(2,j,1)
         r_2(3,nal)=xa(3,j,1)
 10   continue

      call u3b(w,r_1,r_2,nal,1,rms,u,t,ier) !u rotate r_1 to r_2
      print *,'ier=',ier,',rms=',dsqrt(rms)

      do i1=1,nseq1
         coor(1,i1,0)=t(1)+u(1,1)*xa(1,i1,0)+u(1,2)*xa(2,i1,0)+u(1,3)*xa(3,i1,0)
         coor(2,i1,0)=t(2)+u(2,1)*xa(1,i1,0)+u(2,2)*xa(2,i1,0)+u(2,3)*xa(3,i1,0)
         coor(3,i1,0)=t(3)+u(3,1)*xa(1,i1,0)+u(3,2)*xa(2,i1,0)+u(3,3)*xa(3,i1,0)
      enddo

      return
      end

***************************************************************
*     calculate distance between two points from two structures
***************************************************************
      function distp2(xa,moli,molj,i,j)
      include 'pars.h'
      real xa(3,maxr,0:1)
      distp2=(xa(1,i,moli)-xa(1,j,molj))**2
     &      +(xa(2,i,moli)-xa(2,j,molj))**2
     &      +(xa(3,i,moli)-xa(3,j,molj))**2
      return
      end



****************************************************************
*     calculate TM-score and score matrix for DP and LSAP
****************************************************************
      subroutine get_score(invmap,TM,score,score_flag,
     &     mode_flag,quick_flag,mode)
      implicit none
      include 'pars.h'

      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,maxr,0:1)
      common/d0/d0,anseq
      common/contcf/d_col,d_col2
      common/cbvec/cbvec(3,maxr,0:1)  ! normalized vector from ca to cb

      common/id2chain/id_chain1(maxr),id_chain2(maxr)
      common/contlst/cont1,cont2
      common/contnum/ncont1(maxr),ncont2(maxr)  !number of contacts
      integer cont1(maxc,maxr),cont2(maxc,maxr) !contact list
      integer id_chain1,id_chain2,ncont1,ncont2

      integer nseq1,nseq2,nseq,n_al
      integer invmap(maxr),invmap_r(maxr)

      real    xa,d0,anseq,d_col,d_col2,cbvec
      real    xtm1(maxr),ytm1(maxr),ztm1(maxr)
      real    xtm2(maxr),ytm2(maxr),ztm2(maxr)
      real    xtm1r(maxr),ytm1r(maxr),ztm1r(maxr)
      real    xtm2r(maxr),ytm2r(maxr),ztm2r(maxr)
      real    xx,yy,zz,dd
      real    xxb(maxr),yyb(maxr),zzb(maxr)
      real    score(maxr,maxr),product
      real    TM, TMr, Rcomm, diff,diff1,diff2

      real    fcol              !contact overlap factor
      real    d0_input,d0_search
      integer is_ali(maxr,maxr) !original position to aligned position
      integer imap1(maxr),imap2(maxr)
      integer col,ncol,ncolr,Lcomm
      integer score_flag,mode_flag,mode
      integer quick_flag,isearch

ccc   RMSD:
      double precision r_1(3,maxr),r_2(3,maxr),r_1r(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms,drms !armsd is real
      integer ier
      data w /maxr*1.0/
ccc
      integer i,j

      integer seqmat(0:19,0:19)    ! sequence substitution matrix
      integer iseq1(maxr),iseq2(maxr), iaa1,iaa2
      real    sdel,seqsc
      common/seq/iseq1,iseq2,seqmat,sdel


      mode=0
      isearch=2  !extensive search
      if(quick_flag.ge.1) isearch=1

      !if(mode_flag.eq.1) call print_invmap(invmap,nseq2,0)


      !!! deal with single chains here, but could be expanded for multi-chains
      do i=1,nseq1
         id_chain1(i)=1
      enddo
      do j=1,nseq2
         id_chain2(j)=1
      enddo


c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
ccc   for TM-score:
            xtm1(n_al)=xa(1,i,0) !for TM-score
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            imap1(n_al)=i       !aligned position to original index
            imap2(n_al)=j
         endif
      enddo
***   calculate TM-score for the given alignment----------->
      d0_input=d0
      d0_search=d0
      if (d0_search.gt.8.0)d0_search=8.0
      !if (d0_search.lt.3.0)d0_search=3.0

      call TMsearch(d0_input,d0_search,n_al,xtm1,ytm1,ztm1,imap1,
     &     n_al,xtm2,ytm2,ztm2,imap2,TM,Rcomm,Lcomm,isearch,
     &     score_flag,ncol)     !simplified search engine
      TM=TM*n_al/anseq          !TM-score


ccccc invmap_r maps sequence 1 to sequence 2, and use it to
ccccc calculate TM-score again. note that TM-score calculation is
ccccc residue order dependent. The condition block below ensures that
ccccc TM-score is symmetrical during the non-squential alignment search
      if(mode_flag.eq.1)then
         do i=1,nseq1
            invmap_r(i)=-1   !seq1 to seq2, reverse of invmap
         enddo

         do j=1,nseq2
            i=invmap(j)
            if(i.gt.0) invmap_r(i)=j
         enddo

         n_al=0
         do i=1,nseq1
            j=invmap_r(i)         !j aligned to i
            if(j.gt.0)then
               n_al=n_al+1
               xtm1r(n_al)=xa(1,i,0) !for TM-score
               ytm1r(n_al)=xa(2,i,0)
               ztm1r(n_al)=xa(3,i,0)
               xtm2r(n_al)=xa(1,j,1)
               ytm2r(n_al)=xa(2,j,1)
               ztm2r(n_al)=xa(3,j,1)

               r_1r(1,n_al)=xa(1,i,0)
               r_1r(2,n_al)=xa(2,i,0)
               r_1r(3,n_al)=xa(3,i,0)

               imap1(n_al)=i    !aligned position to original index
               imap2(n_al)=j
            endif
         enddo

         call TMsearch(d0_input,d0_search,n_al,xtm1r,ytm1r,ztm1r,
     &        imap1,n_al,xtm2r,ytm2r,ztm2r,imap2,TMr,Rcomm,Lcomm,
     &        isearch,score_flag,ncolr)  !simplified search engine
         TMr=TMr*n_al/anseq       !TM-score

         !write(*,'(2(a,F8.5))')'TM=',TM,' TMr=',TMr

         !!! ensure the alignment is same as that from two structures in revsersed order
         diff=TMr-TM
         diff1=abs(xa(1,1,0)-xa(1,2,0)) ! arbitrarily break symmetry
         diff2=abs(xa(1,1,1)-xa(1,2,1))
         if(diff.gt.1e-6.or.(abs(diff).le.1e-6.and.diff1.gt.diff2))then
            TM=TMr
            ncol=ncolr
            mode=1
            do i=1,n_al
               xtm1(i)=xtm1r(i)
               ytm1(i)=ytm1r(i)
               ztm1(i)=ztm1r(i)

               r_1(1,i)=r_1r(1,i)
               r_1(2,i)=r_1r(2,i)
               r_1(3,i)=r_1r(3,i)
            enddo
         endif
      endif




ccc   is_ali(i,j) indicate whether i, j are in alignment dis<d_col
      do i=1,nseq1
         do j=1,nseq2
            is_ali(i,j) = -1
         enddo
      enddo


***   calculate score matrix score(i,j)------------------>
      do i=1,n_al
         r_2(1,i)=xtm1(i)
         r_2(2,i)=ytm1(i)
         r_2(3,i)=ztm1(i)
      enddo
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      !call print_TransMatrix(t,u)

      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)

         xxb(i)=u(1,1)*cbvec(1,i,0)+u(1,2)*cbvec(2,i,0)+u(1,3)*cbvec(3,i,0)
         yyb(i)=u(2,1)*cbvec(1,i,0)+u(2,2)*cbvec(2,i,0)+u(2,3)*cbvec(3,i,0)
         zzb(i)=u(3,1)*cbvec(1,i,0)+u(3,2)*cbvec(2,i,0)+u(3,3)*cbvec(3,i,0)
         do j=1,nseq2
            if(id_chain1(i).eq.id_chain2(j))then
               dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
               score(i,j)=1/(1+dd/d0**2)

               if(dd.le.d_col2)then
                  is_ali(i,j)=1
               endif
            else
               score(i,j)=ninf
            endif
         enddo
      enddo

ccc   calculate contact overlap and get IS score
      ncol=0  ! total number of contact overlaps
      do i=1,nseq1
         do j=1,nseq2
            if(score_flag.eq.1.and.id_chain1(i).eq.id_chain2(j))then
               call get_fcol_4s(i,j,is_ali,fcol,col,
     &              cont1,cont2,ncont1,ncont2)
               ncol=ncol+col
               score(i,j)=score(i,j)*(fcol+0.01)
            else if(score_flag.eq.2)then
               iaa1=iseq1(i)
               iaa2=iseq2(j)
               if(iaa1.ge.0.and.iaa2.ge.0)then
                  seqsc=seqmat(iaa1,iaa2)
                  if(seqsc.gt.0)then
                     score(i,j)=score(i,j)*seqsc
                  else
                     score(i,j)=score(i,j)*sdel
                  endif
                  if(iaa1.gt.0.and.iaa2.gt.0)then
                     product=xxb(i)*cbvec(1,j,1)
                     product=product+yyb(i)*cbvec(2,j,1)
                     product=product+zzb(i)*cbvec(3,j,1)
                  else if(iaa1.eq.0.and.iaa2.eq.0)then
                     product=1   !both are glycine
                  else
                     product=0.272 ! mean cos(theta) from random alignments
                  endif
                  product=product+0.5
                  if(product.gt.1)then
                     product=1
                  else if(product.lt.0.1)then
                     product=0.1
                  endif
                  score(i,j)=score(i,j)*product
               else
                  score(i,j)=score(i,j)*sdel
               endif
            endif
            !print *,i,j,score(i,j)
         enddo
      enddo
      ncol=ncol/2  !each contact counted twice

      !print *,'get_score is done'
c^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


****************************************************************
*     with invmap(i) calculate score(i,j) using RMSD rotation
****************************************************************
      subroutine get_score1(invmap)
      include 'pars.h'
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,maxr,0:1)
      common/dpc/score(maxr,maxr),gap_open
      common/d0/d0,anseq
      common/d0min/d0_min
      dimension xtm1(maxr),ytm1(maxr),ztm1(maxr)
      dimension xtm2(maxr),ytm2(maxr),ztm2(maxr)

      integer invmap(maxr)


ccc   RMSD:
      double precision r_1(3,maxr),r_2(3,maxr),r_3(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /maxr*1.0/
ccc


c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
         endif
      enddo
***   calculate score matrix score(i,j)------------------>
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      d01=d0+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01*d01
      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq2
            dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
            score(i,j)=1/(1+dd/d02)
         enddo
      enddo

c^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end






********************************************************************
*     Dynamic programming for alignment.
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
*
*     Please note this subroutine is not a correct implementation of
*     the N-W dynamic programming because the score tracks back only
*     one layer of the matrix. This code was exploited in TM-align
*     because it is about 1.5 times faster than a complete N-W code
*     and does not influence much the final structure alignment result.
********************************************************************
      SUBROUTINE DP(score,gap_open,NSEQ1,NSEQ2,invmap)
      include 'pars.h'
      LOGICAL*1 DIR
      real score(maxr,maxr),gap_open
      integer invmap(maxr)
      dimension DIR(0:maxr,0:maxr),VAL(0:maxr,0:maxr)
      common/speed/align_len1
      integer align_len1

      REAL H,V
      align_len1=0
***   initialize the matrix:
      val(0,0)=0
      do i=1,nseq1
        dir(i,0)=.false.
        val(i,0)=0
      enddo
      do j=1,nseq2
        dir(0,j)=.false.
        val(0,j)=0
        invmap(j)=-1
      enddo

***   decide matrix and path:
      DO j=1,NSEQ2
        DO i=1,NSEQ1
          D=VAL(i-1,j-1)+SCORE(i,j)
          H=VAL(i-1,j)
          if(DIR(i-1,j))H=H+GAP_OPEN
          V=VAL(i,j-1)
          if(DIR(i,j-1))V=V+GAP_OPEN

          IF((D.GE.H).AND.(D.GE.V)) THEN
            DIR(I,J)=.true.
            VAL(i,j)=D
          ELSE
            DIR(I,J)=.false.
            if(V.GE.H)then
              val(i,j)=v
            else
              val(i,j)=h
            end if
          ENDIF
        ENDDO
      ENDDO

***   extract the alignment:
      i=NSEQ1
      j=NSEQ2
      DO WHILE((i.GT.0).AND.(j.GT.0))
        IF(DIR(i,j))THEN
          invmap(j)=i
          align_len1=align_len1+1
          i=i-1
          j=j-1
        ELSE
          H=VAL(i-1,j)
          if(DIR(i-1,j))H=H+GAP_OPEN
          V=VAL(i,j-1)
          if(DIR(i,j-1))V=V+GAP_OPEN
          IF(V.GE.H) THEN
            j=j-1
          ELSE
            i=i-1
          ENDIF
        ENDIF
      ENDDO

c^^^^^^^^^^^^^^^Dynamical programming done ^^^^^^^^^^^^^^^^^^^
      return
      END



****************************************************************
*     remove Ca pairs with distance larger than d8 from invmap
****************************************************************
      subroutine remove_distant_pairs(invmap,nseq1,nseq2,xa)
      implicit none
      include 'pars.h'

      integer nseq1,nseq2,nseq
      integer invmap(maxr),invmap1(maxr)

      common/d8/d8
      real    d8

      real    xtm1(maxr),ytm1(maxr),ztm1(maxr)
      real    xtm2(maxr),ytm2(maxr),ztm2(maxr)
      real    xx,yy,zz
      real    XA(3,maxr,0:1)
      real    dd(maxr),dcut

ccc   RMSD:
      double precision r_1(3,maxr),r_2(3,maxr),r_3(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms,drms !armsd is real
      integer ier
      data w /maxr*1.0/
ccc
      integer i,j,n_al

      dcut=d8
c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,nseq2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
         endif
      enddo
***   calculate score matrix score(i,j)------------------>
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2

      if(n_al.le.3)return

      do 10 j=1,nseq2
         i=invmap(j)
         if(i.le.0) goto 10
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)

         dd(j)=sqrt((xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2)
 10   continue

 20   n_al=0
      do 30 j=1,nseq2
         i=invmap(j)
         if(i.gt.0.and.dd(j).le.dcut) then
            invmap1(j) = i
            n_al=n_al+1
         else
            invmap1(j) = -1
         endif
 30   continue
      if(n_al.lt.3)then
         dcut=dcut+1
         goto 20
      endif

      do j=1,nseq2
         invmap(j)=invmap1(j)
      enddo

      return
      end



********************************************************************
*     Sequence order dependent matching
*
*     Use an efficient algorithm for solving a linear sum assignment
*     problem at time complexity of O(n^3).
*
*     Mu Gao
********************************************************************
      SUBROUTINE solvLSAP(score,invmap,nseq1,nseq2,n,z)
      implicit none
      include 'pars.h'

      integer invmap(maxr)
      integer nseq1, nseq2
      integer n
      integer spalte(n)
      integer i,j

      real c(n,n), cost, z
      real eps,sup
      real score(maxr,maxr)

      character*50 filename

      eps = 1E-10
      sup = 1E+9


      ! cost matrix for linear sum assignment problem
      do i = 1, n
        do j = 1, n
           c(j,i) = 2
        enddo
      enddo

      do i = 1, n
        do j = 1, n
           if(i.le.nseq1.and.j.le.nseq2)then
              c(i,j) = 2 - score(i,j)
           endif
        enddo
      enddo



      call lsapr ( n, c, z, spalte, sup, eps )

      do j = 1, nseq2
         i=spalte(j)
         if(i.le.nseq1)then
            invmap(j) = i
         else
            invmap(j) = -1
         endif
      end do
      !write(*,'(a,f8.4)')'optimal cost z=',z,

      return
      END
*****************  End of solvLSAP subroutein **************************




      subroutine lsapr ( n, cost, z, spalte, sup, eps )

c*********************************************************************72
c
cc LSAPR solves the linear sum assignment problem with real data.
c
c  Modified:
c
c    22 November 2007
c
c  Author:
c
c    Ulrich Derigs
c
c  Modified by Mu Gao
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the dimension of the cost matrix.
c
c    Input, real SUP, a large machine number.
c
c    Input, real C(N,N), the cost matrix.
c
c    Output, real Z, the optimal value.
c
c    Output, integer SPALTE(N), the optimal assignment.
c
c    Input, real EPS, machine accuracy.
c

      implicit none

      integer n

      real cost(n,n)
      real cc
      real d    !Mu: this should be real!
      real dminus(n)
      real dplus(n)
      real eps,sup
      real ui
      real vgl
      real vj
      real ys(n)
      real yt(n)
      real z

      integer i,j,u,w
      integer ind,index,j0
      integer zeile(n)
      integer vor(n)
      integer spalte(n)

      logical label(n)

c
c  Construct an initial partial assignment.
c
      do i = 1, n
        zeile(i) = 0
        spalte(i) = 0
        vor(i) = 0
        ys(i) = 0.0E+00
        yt(i) = 0.0E+00
      end do

      do i = 1, n
        do j = 1, n
          cc = cost(j,i)

          if ( j .ne. 1 ) then
            if ( ( cc - ui ) .ge. eps ) then
              go to 3
            end if
          end if

          ui = cc
          j0 = j

3         continue
        end do

        ys(i) = ui

        if ( zeile(j0) .eq. 0 ) then
          zeile(j0) = i
          spalte(i) = j0
        end if

      end do

      do j = 1, n
        yt(j) = 0
        if ( zeile(j) .eq. 0 ) then
          yt(j) = sup
        end if
      end do


      do i = 1, n
        ui = ys(i)

        do j = 1, n
          vj = yt(j)
          if ( eps .lt. vj ) then

            cc = cost(j,i) - ui
            if ( cc + eps .lt. vj ) then
              yt(j) = cc
              vor(j) = i
            end if

          end if
        end do
      end do

      do j = 1, n
        i = vor(j)
        if ( i .ne. 0 ) then
          if ( spalte(i) .eq. 0 ) then
            spalte(i) = j
            zeile(j) = i
          end if
        end if
      end do

      do i = 1, n
        if ( spalte(i) .eq. 0 ) then
          ui = ys(i)

          do j = 1, n
            if ( zeile(j) .eq. 0 ) then
              cc = cost(j,i)
              if ( cc - ui - yt(j) .le. -eps ) then
              !if ( ( cc - ui - yt(j) + eps ) .le. 0.0E+00 ) then
                spalte(i) = j
                zeile(j) = i
              end if
            end if
          end do

        end if
      end do
c
c  Construct the optimal solution.
c
      do 1000 u = 1, n

        if ( spalte(u) .gt. 0 ) goto 1000

c
c  Shortest path computation.
c

        do i = 1, n
          vor(i) = u
          label(i) = .false.
          dplus(i) = sup
          dminus(i) = cost(i,u) - ys(u) - yt(i)
        end do

        dplus(u) = 0.0E+00

105     continue

        d = sup
        index = 0

        do i = 1, n
          if ( .not. label(i) ) then
            if ( dminus(i) + eps .lt. d ) then
              d = dminus(i)
              index = i
            end if
          end if
        end do

        if ( index .eq. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LSAPR - Fatal error!'
          write ( *, '(a)' ) '  No unlabeled node with DMINUS < D.'
          stop
        end if

        if ( zeile(index) .le. 0 ) then
          go to 400
        end if

        label(index) = .true.
        w = zeile(index)
        dplus(w) = d

        do i = 1, n
          if ( .not. label(i) ) then
            vgl = d + cost(i,w) - ys(w) - yt(i)
            if ( vgl + eps .lt. dminus(i) ) then
              dminus(i) = vgl
              vor(i) = w
            end if
          end if
        end do

        go to 105
c
c  Augmentation.
c
400     continue

        w = vor(index)
        zeile(index) = w
        ind = spalte(w)
        spalte(w) = index

        if ( w .ne. u ) then
          index = ind
          go to 400
        end if
c
c  Transformation.
c
!500     continue

        do i = 1, n
          if ( dplus(i) .ne. sup ) then
            ys(i) = ys(i) + d - dplus(i)
          end if

          if ( dminus(i) + eps .lt. d ) then
            yt(i) = yt(i) + dminus(i) - d
          end if
        end do

1000  continue
c
c  Computation of the optimal value.
c
      z = 0.0E+00
      do i = 1, n
        j = spalte(i)
        z = z + cost(j,i)
      end do

      return
      end
*****************  End of LSAPR subroutein **************************







ccccc================================================ccccc
ccccc print the mapping array for debugging purpose  ccccc
ccccc================================================ccccc
      subroutine print_invmap(invmap,len,id)
      implicit none
      include 'pars.h'

      common/initial4/mm1(maxr),mm2(maxr)
      integer invmap(maxr),len,mm1,mm2
      integer i,j,prev_i,prev_j,id

      print *,'printing invmap...'

      prev_i=-1
      prev_j=-1
      do j=1,len
         i=invmap(j)
         if(i.gt.0)then
            if(i-prev_i.ne.1.or.j-prev_j.ne.1)then
               write(*,'(A)') '--'
            endif
            !write(*,'(2I6)') i, j
            if(id.eq.1)then
               write(*,'(2I6,A,2I6)') mm2(i), mm1(j), ', ind: ', i,j
            else
               write(*,'(2I6,A,2I6)') mm1(i), mm2(j), ', ind: ', i,j
            endif
         endif
         prev_i=i
         prev_j=j
      enddo
      write(*,'(/)')

      return
      end


ccccc================================================ccccc
ccccc   reverse invmap from seq1 -> seq2 to seq2 -> seq1
ccccc================================================ccccc
      subroutine rev_invmap(invmap_r,invmap,nseq1,nseq2)
      implicit none
      include 'pars.h'

      integer invmap(maxr),invmap_r(maxr)
      integer i,j,nseq1,nseq2

      do i=1,nseq2
         invmap(i)=-1
      enddo

      do i=1,nseq1
         j=invmap_r(i)
         if(j.gt.0) invmap(j)=i            !seq2 to seq1
      enddo

      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Calculate Ca-Ca contact using a distance cutoff                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calCaCont(coor,id,nseq,ncont,cont,cutoff,dcont)
      implicit none
      include 'pars.h'

      real coor(3,maxr,0:1)  ! ca cooridnates
      real dist, dist2       ! distance squared
      real cutoff,cutoff2
      real dis(maxc,maxr)   ! distance matrix
      real sd(maxc),sc(maxc)
      real dcont(maxr,maxr)

      integer ncont(maxr),cont(maxc,maxr)
      integer id             ! id of the structure
      integer nseq           ! length of structure
      integer i,j,k

      cutoff2=cutoff*cutoff

****** Initializing contact arrays ********
      do i=1,maxr
         ncont(i)=0
         do j=1,maxc
            cont(j,i)=-1
         enddo
      enddo

****** Calculating Ca-Ca contacts ********
      do 10 i=1,nseq
         do 20 j=i+1,nseq
            dist2 = (coor(1,i,id)-coor(1,j,id))**2 +
     &              (coor(2,i,id)-coor(2,j,id))**2 +
     &              (coor(3,i,id)-coor(3,j,id))**2
            dist=sqrt(dist2)
            dcont(i,j)=dist
            dcont(j,i)=dist
            if( dist .lt. cutoff ) then
               if(ncont(i).eq.maxc.or.ncont(j).eq.maxc) goto 20
               ncont(i) = ncont(i) + 1
               ncont(j) = ncont(j) + 1
               cont(ncont(i),i) = j
               cont(ncont(j),j) = i
               dis(ncont(i),i) = dist
               dis(ncont(j),j) = dist
            endif
 20      enddo
 10   enddo

      do i=1,nseq
         do j=1,ncont(i)
            k=cont(j,i)
            sd(j)=dis(j,i)
            sc(j)=k
            !write(*,'(A,3I6,F8.3)')'Dist',i,j,k,dis(j,i)
         enddo

         if(ncont(i).gt.1) call ssort( sd, sc, ncont(i), 2 )

         do j=1,ncont(i)
            dis(j,i)=sd(j)
            cont(j,i)=nint(sc(j))
            !write(*,'(A,3I6,F8.3)')'Dist',i,j,cont(j,i),dis(j,i)
         enddo

      enddo



      return
      end


      SUBROUTINE SSORT (X, Y, N, KFLAG)
C***BEGIN PROLOGUE  SSORT
C***PURPOSE  Sort an array and optionally make the same interchanges in
C            an auxiliary array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2B
C***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   SSORT sorts array X and optionally makes the same interchanges in
C   array Y.  The array X may be sorted in increasing order or
C   decreasing order.  A slightly modified quicksort algorithm is used.
C
C   Description of Parameters
C      X - array of values to be sorted   (usually abscissas)
C      Y - array to be (optionally) carried along
C      N - number of values in array X to be sorted
C      KFLAG - control parameter
C            =  2  means sort X in increasing order and carry Y along.
C            =  1  means sort X in increasing order (ignoring Y)
C            = -1  means sort X in decreasing order (ignoring Y)
C            = -2  means sort X in decreasing order and carry Y along.
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***REVISION HISTORY  (YYMMDD)
C   761101  DATE WRITTEN
C   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced statement labels.  (WRB)
C   891024  Changed category.  (WRB)
C   891024  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   901012  Declared all variables; changed X,Y to SX,SY. (M. McClain)
C   920501  Reformatted the REFERENCES section.  (DWL, WRB)
C   920519  Clarified error messages.  (DWL)
C   920801  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
C***END PROLOGUE  SSORT
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      REAL X(*), Y(*)
C     .. Local Scalars ..
      REAL R, T, TT, TTY, TY
      INTEGER I, IJ, J, K, KK, L, M, NN
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. External Subroutines ..
C     None
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***FIRST EXECUTABLE STATEMENT  SSORT
      NN = N
      IF (NN .LT. 1) THEN
         PRINT *,
     +      'The number of values to be sorted is not positive.'
         RETURN
      ENDIF
C
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         PRINT *,
     +      'The sort control parameter, K, is not 2, 1, -1, or -2.'
         RETURN
      ENDIF
C
C     Alter array X to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            X(I) = -X(I)
   10    CONTINUE
      ENDIF
C
      IF (KK .EQ. 2) GO TO 100
C
C     Sort X only
C
      M = 1
      I = 1
      J = NN
      R = 0.375E0
C
   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
   30 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = X(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than than T, interchange with T
C
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
   40 L = L-1
      IF (X(L) .GT. T) GO TO 40
C
C     Find an element in the first half of the array which is greater
C     than T
C
   50 K = K+1
      IF (X(K) .LT. T) GO TO 50
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         GO TO 40
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70
C
C     Begin again on another portion of the unsorted array
C
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
C
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = X(I+1)
      IF (X(I) .LE. T) GO TO 80
      K = I
C
   90 X(K+1) = X(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 90
      X(K+1) = T
      GO TO 80
C
C     Sort X and carry Y along
C
  100 M = 1
      I = 1
      J = NN
      R = 0.375E0
C
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
  120 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = X(IJ)
      TY = Y(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
         Y(IJ) = Y(I)
         Y(I) = TY
         TY = Y(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than T, interchange with T
C
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
         Y(IJ) = Y(J)
         Y(J) = TY
         TY = Y(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
            Y(IJ) = Y(I)
            Y(I) = TY
            TY = Y(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
  130 L = L-1
      IF (X(L) .GT. T) GO TO 130
C
C     Find an element in the first half of the array which is greater
C     than T
C
  140 K = K+1
      IF (X(K) .LT. T) GO TO 140
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         TTY = Y(L)
         Y(L) = Y(K)
         Y(K) = TTY
         GO TO 130
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
C
C     Begin again on another portion of the unsorted array
C
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
C
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = X(I+1)
      TY = Y(I+1)
      IF (X(I) .LE. T) GO TO 170
      K = I
C
  180 X(K+1) = X(K)
      Y(K+1) = Y(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 180
      X(K+1) = T
      Y(K+1) = TY
      GO TO 170
C
C     Clean up
C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            X(I) = -X(I)
  200    CONTINUE
      ENDIF
      RETURN
      END


ccccc====================================================cccccc
      subroutine init_ali_map( is_ali, nseq1, nseq2 )
      include 'pars.h'

      integer nseq1,nseq2
      integer is_ali(maxr, maxr)   !alignment map
      integer i,j

ccc   is_ali(i,j) indicate whether i, j are in alignment dis<d_col
      do i=1,nseq1
         do j=1,nseq2
            is_ali(i,j) = -1
         enddo
      enddo

      return
      end




ccccc===================================================ccccc
ccccc  print transformation matrix                      ccccc
ccccc===================================================ccccc
      subroutine print_TransMatrix( t, u )

      real*8 u(3,3)
      real*8 t(3)

      ! overwrite old matrix
      do i=1,3
         write(*,'(I5,3F6.2f)')nint(t(i)),(u(i,j),j=1,3)
      enddo
      write(*,*)

      return
      end


ccccc===================================================ccccc
ccccc  Print the help message                           ccccc
ccccc===================================================ccccc
      subroutine printHelp()

         write(*,*)
         write(*,'(A)')'Usage: apoc <options> pdbfile1 pdbfile2 '
         write(*,*)
         write(*,'(A)')'Options:'

         write(*,'(A)')'Input options'
         write(*,'(A,/,A)')' -fa <num>', '    Global structure'//
     &     ' alignment: 1 - enable (default), 0 - disable.'

         write(*,'(A,/,A)')' -lt <file>','    Provide a'//
     &     ' list of templates to compare in a file.'
         write(*,'(A,/,A)')' -lq <file>','    Provide a'//
     &     ' list of queries (targets) to compare in a file.'

         write(*,'(A,/,A)')' -pt <str1,str2,...>','    Names of pockets'//
     &     ' in the first (template) structure for comparison.'
         write(*,'(A,/,A)')' -pq <str1,str2,...>','    Names of pockets'//
     &     ' in the second (query) structure for comparison.'
         write(*,'(A,/,A)')' -block <file>','    Load a block of '//
     &     'concatenated pdb files.'

         write(*,'(A,/,A)')' -pvol <num>','    Minimal pocket volume'//
     &     ' in grid points. Default 100'
         write(*,'(A,/,A)')' -plen <num>','    Minimal number of'//
     &     ' pocket residues. Default 10'

         write(*,'(/,A)')'Alignment options'
         write(*,'(A,/,A)')' -sod','    Restrict to sequence-order-dependent'//
     &        ' alignment. Default no restriction.'
         write(*,'(A,/,A)')' -v','    Alignment printout: 0 - none,'//
     &     ' 1 - concise, 2 - detailed (default).'

         write(*,'(/,A)')'Scoring options'
         write(*,'(A,/,A)')' -m <str>',  '    Similarity scoring metric: '//
     &        ' tm (TM-score), ps (PS-score, default).'

         write(*,'(A,/,A)')' -L <num>','    Normalize the score with a fixed '//
     &        'length specified by num.'

         write(*,'(A,/,A)')' -a','    Normalize the score by the average'//
     &        ' size of two structures.'
         write(*,'(A,/,A)')' -b','    Normalize the score by the minimum'//
     &        ' size of two structures.'
         write(*,'(A,/,A)')' -c','    Normalize the score by the maximum'//
     &        ' size of two structures.'

         write(*,*)

      return
      end

ccccc==========================================================cccc


      subroutine getBasename(fullname,basename)
      implicit none
      character fullname*(*),basename*(*)
      integer i,j,l,n_sta,n_end

      l = len(fullname)
      n_sta=0
      n_end=l
      do i=1,l
         if(fullname(i:i).eq.'/')then
            n_sta=i
         else if(fullname(i:i).ne.' ')then
            n_end=i
         endif
      enddo

      l = n_end - n_sta
      do i=1,l
         j=n_sta+i
         basename(i:i)=fullname(j:j)
      enddo

      ! maximum length of basename is 30 characters
      do i=l+1,30
         basename(i:i)=' '
      enddo

      return
      end


ccccc==========================================================cccc
ccccc read specified pockets for comparison
ccccc==========================================================cccc
      subroutine getPKind(fullname,pksel,npksel)
      use Constants
      implicit none
      character fullname*(*)
      character*30 pksel(maxp)
      integer npksel, i,j,l,n_sta,n_end

      l = len(fullname)
      n_sta=1
      n_end=0
      npksel=0
      do i=1,l
         if(fullname(i:i).eq.',')then
            n_end=i-1
            npksel=npksel+1
            read(fullname(n_sta:n_end),*)pksel(npksel)
            n_sta=i+1
         else if(fullname(i:i).eq.' ')then
            n_end=i-1
            npksel=npksel+1
            read(fullname(n_sta:n_end),*)pksel(npksel)
            exit
         endif
      enddo

      return
      end




ccccc==========================================================cccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Read a single PDB file                                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readPDBFileNew(pdbfile,pdb)
      use PDBType
      use Constants
      implicit none

      character*500 pdbfile
      character*100 line, pdblines(maxs)
      character*3 head
      integer nline,nres,npkt,flag

      type(pdbstruct) :: pdb

      pdb%filename = pdbfile

      nline=0
      nres=0
      npkt=0
      flag=1
      open(unit=8,file=pdbfile,status='old')
      do while (.true.)
         read(8,'(A100)',end=100)line
         head=line(1:3)
         if(head.eq.'ATO'.or.head.eq.'HET'.or.head.eq.'TER'.or.head.eq.'PKT')then
            nline=nline+1
            pdblines(nline)=line
         endif
         if(line(1:4).eq.'ATOM'.and.(line(13:16).eq.'CA  '.or.
     &      line(13:16).eq.' CA '.or.line(13:16).eq.'  CA').and.flag.eq.1) nres=nres+1
         if(head.eq.'PKT')then
            npkt=npkt+1
            flag=0
         endif
      enddo
 100  continue
      close(8)

      call readPDBFromLines(pdb,pdblines,nline,nres,npkt)

      return
      end

ccccc==========================================================cccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Read concatenated PDB records in one file                        c
c     Each PDB record starts with "PDBSTA name" and ends with "PDBEND" c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  readPDBBlock(block_file,struct_num,struct)
      use Constants
      use PDBType
      implicit none

      type(pdbstruct) :: struct(maxt)

      integer struct_num, nline, nres, npkt, flag
      character*500 block_file,pdb_file_name
      character*100 line, pdblines(maxs)
      character*6 head

      open(unit=8,file=block_file,status='old')

      do while (.true.)
         read(8,'(A100)',end=100)line
         if(line(1:6).eq.'PDBSTA')then
            struct_num=struct_num+1
            nline=0
            nres=0
            npkt=0
            flag=1
            read(line,'(A6,A)')head,pdb_file_name
            struct(struct_num)%filename=adjustl(pdb_file_name)
         else if(line(1:6).eq.'PDBEND')then
            call readPDBFromLines(struct(struct_num),pdblines,nline,nres,npkt)
         else
            nline=nline+1
            pdblines(nline)=line
            if(line(1:4).eq.'ATOM'.and.(line(13:16).eq.'CA  '.or.
     &         line(13:16).eq.' CA '.or.line(13:16).eq.'  CA').and.flag.eq.1)nres=nres+1
            if(line(1:3).eq.'PKT') then
               npkt=npkt+1
               flag=0
            endif
         endif
      enddo
 100  continue
      close(8)


      return
      end
ccccc==========================================================cccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Read a PDB record.   Please check the PDB format at:             c
c     http://www.wwpdb.org/docs.html#format                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readPDBFromLines( pdb, pdblines, num_lines, nres, npkt )
      use PDBType
      implicit none

      character*100 line, pdblines(maxs)
      character*1 altloc,chain,chtmp,ic
      character*3 resnam,chtmp3
      character*5 atmnam
      character*4 resid
      character*6 head
      character*7 res
      character*7,pointer :: rescb(:),resca(:)
      character*30 pname

      real,pointer ::  coor_tmp(:,:)
      real   x,y,z,occ,temp,r1,r2,r3,r

      integer  nres,npkt
      integer  nchain, ichterm(0:maxk)
      integer  length,length2     !total number of residues
      integer  iatom,ires
      integer  i,j,k,kk,nc,counter
      integer  nline,num_lines
      integer  num_pkt,index,pkt_flag,num_pkt_res
      integer  len,vol

      character*3 aa(-1:19)
      data aa/ 'UNK','GLY','ALA','SER','CYS',
     &     'VAL','THR','ILE','PRO','MET',
     &     'ASP','ASN','LEU','LYS','GLU',
     &     'GLN','ARG','HIS','PHE','TYR',
     &     'TRP'/

      character*1 slc(-1:19)
      data slc/'X','G','A','S','C',
     &     'V','T','I','P','M',
     &     'D','N','L','K','E',
     &     'Q','R','H','F','Y',
     &     'W'/

      type(pdbstruct) :: pdb

      if(nres.lt.0)then
         print '(A)','Error: no protein residue found in your structure!'
         stop
      endif

      if(nres.gt.maxr)then
         print '(A,I)','Warning: total number of protein residues is more than the maximum',maxr
      endif

      allocate(pdb%resname(nres))
      allocate(pdb%seq(0:nres))
      allocate(pdb%icode(nres))
      allocate(pdb%chain(nres))
      allocate(pdb%resid(nres))
      allocate(pdb%sec_struct(nres))
      allocate(pdb%iseq(nres))
      allocate(pdb%coor_ca(3,nres))
      allocate(pdb%coor_cb(3,nres))
      allocate(pdb%cbvec(3,nres))
      allocate(pdb%pk_res(npkt,nres))
      allocate(pdb%pk_len(npkt))
      allocate(pdb%pk_vol(npkt))
      allocate(pdb%pk_name(npkt))

      allocate(coor_tmp(3,nres))
      allocate(rescb(nres))
      allocate(resca(nres))

      i=0
      k=0
      kk=0
      nchain=0
      chain=' '
      ichterm(0)=0
      num_pkt=0
      pkt_flag=0
      nline=0
      do nline=1,num_lines
         line=pdblines(nline)

         if(line(1:3).eq.'PKT')then
            pkt_flag=1
            num_pkt=num_pkt+1
            pdb%pk_len(num_pkt)=0
            read(line,*)head,len,vol,pname
            !print '(A,3I8,1X,A)','found pocket ',num_pkt,len,vol,pname
            pdb%pk_vol(num_pkt)=vol
            pdb%pk_name(num_pkt)=pname
         else if(line(1:3).eq.'TER'.and.pkt_flag.eq.0)then
            nchain=nchain+1     ! count number of chains
            ichterm(nchain)=i   ! index of the last residue of a chain
         endif

         !!! read full length residues
         if(line(1:4).eq.'ATOM'.and.pkt_flag.eq.0)then
            read(line,9000) head,iatom,chtmp,atmnam,altloc,resnam,
     &           chtmp,chain,resid,ic,chtmp3,x,y,z,occ,temp

            res=chain//'.'//resid//ic

            if(altloc.ne.' '.and.altloc.ne.'A'.and.altloc.ne.'1')
     &           goto 90

            if(atmnam.eq.'CA  '.or.atmnam.eq.' CA '
     &           .or.atmnam.eq.'  CA') then
               i=i+1
               pdb%coor_ca(1,i)=x
               pdb%coor_ca(2,i)=y
               pdb%coor_ca(3,i)=z
               pdb%resname(i)=resnam
               read(resid,'(I4)')pdb%resid(i)
               pdb%icode(i)=ic
               pdb%chain(i)=chain
               resca(i)=res

               do j=-1,19
                  if(resnam.eq.aa(j))then
                     pdb%seq(i)=slc(j)
                     pdb%iseq(i)=j
                     goto 21
                  endif
               enddo
               pdb%seq(i)=slc(-1)
               pdb%iseq(i)=-1
 21            continue

            !print "(I5,A2,A5,A6,I4)",i,pdb%seq(i),atmnam,pdb%resname(i),pdb%resid(i)
            else if(atmnam.eq.'CB  '.or.atmnam.eq.' CB '
     &           .or.atmnam.eq.'  CB') then
               k=k+1
               rescb(k)=res
               coor_tmp(1,k)=x
               coor_tmp(2,k)=y
               coor_tmp(3,k)=z
            endif

            if(i.ge.maxr)then
               print *, "Warning: total number of residues exceeds ",maxr
               goto 100
            endif

         !!! read pocket residues
         else if(line(1:4).eq.'ATOM'.and.pkt_flag.eq.1)then
            read(line,9000) head,iatom,chtmp,atmnam,altloc,resnam,
     &           chtmp,chain,resid,ic,chtmp3,x,y,z,occ,temp

            res=chain//'.'//resid//ic

            if(altloc.ne.' '.and.altloc.ne.'A'.and.altloc.ne.'1')
     &           goto 90
            if(atmnam.eq.'CA  '.or.atmnam.eq.' CA '
     &           .or.atmnam.eq.'  CA') then
               kk=pdb%pk_len(num_pkt)+1
               do j=1,i
                  if(resca(j).eq.res)then
                     pdb%pk_res(num_pkt,kk)=j
                     pdb%pk_len(num_pkt)=kk
                     !print '(I5,A10,I4,A2,I4)',j,res,kk,pdb%chain(j),pdb%resid(j)
                  endif
               enddo
            endif
         endif

 90      continue
      enddo
 100  continue
 9000 format(A6,I5,A1,A4,A1,A3,2A1,A4,A1,A3,3F8.3,2F6.2)


      !print '(4(A,I6))','nres=',nres,', npkt=',npkt, ', i=',i,', num_pkt=',num_pkt
      pdb%num_residue=i
      pdb%num_pocket=num_pkt

      call calcSecStruct( pdb )


      length2=k
      k=1
      do i=1,pdb%num_residue
!         print '(I5,A4,A10,3F8.3)',i,pdb%seq(i),resca(i),pdb%coor_ca(1,i),pdb%coor_ca(2,i),
!     &        pdb%coor_ca(3,i)
         x=pdb%coor_ca(1,i)
         y=pdb%coor_ca(2,i)
         z=pdb%coor_ca(3,i)
         do j=k,length2
            if(resca(i).eq.rescb(j))then
               x=coor_tmp(1,j)
               y=coor_tmp(2,j)
               z=coor_tmp(3,j)
               k=j
               goto 105
            endif
         enddo
 105     continue
         pdb%coor_cb(1,i)=x
         pdb%coor_cb(2,i)=y
         pdb%coor_cb(3,i)=z

         !!Calculate a normalized vector from Ca to Cb
         r1=pdb%coor_cb(1,i)-pdb%coor_ca(1,i)
         r2=pdb%coor_cb(2,i)-pdb%coor_ca(2,i)
         r3=pdb%coor_cb(3,i)-pdb%coor_ca(3,i)
         r=sqrt(r1*r1+r2*r2+r3*r3)
         if(r.gt.0.1)then
            pdb%cbvec(1,i)=r1/r
            pdb%cbvec(2,i)=r2/r
            pdb%cbvec(3,i)=r3/r
         else
            pdb%cbvec(1,i)=1
            pdb%cbvec(2,i)=0
            pdb%cbvec(3,i)=0
         endif
!         print "(I5,A4,A10,6F8.3,I4)",i,pdb%seq(i),resca(i),pdb%coor_ca(1,i),pdb%coor_ca(2,i),
!     &        pdb%coor_ca(3,i),pdb%coor_cb(1,i),pdb%coor_cb(2,i),pdb%coor_cb(3,i),pdb%sec_struct(i)
      enddo

      deallocate(coor_tmp)
      deallocate(rescb)
      deallocate(resca)

      return
      end

ccccc==========================================================cccc







cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Read an alignment from a file                                    c
c     The format of file is                                            c
c     Index Ch1 Resid1 - AA1 Ch2 Resid2                                c
c     e.g.                                                             c
c
cPREALN 13
c    1 A   385 R - A   407 R
c    2 A   387 S - A   409 S
c    3 A   388 G - A   410 G
c    4 A   390 D - A   412 A
c    5 A   393 V - A   415 I
c    6 A   394 K - A   416 Q
c    7 A   420 S - A   444 S
c    8 A   423 K - A   447 K
c    9 A   424 D - A   448 L
c   10 A   427 R - A   451 R
c   11 A   499 V - A   524 V
c   12 A   502 P - A   527 P
c   13 A   503 T - A   528 N
cEND
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine readAlnFile(matchfile,prealn,nprealn,prealn_len,
     &     chain1,mm1,icode1,seq1,nseq1,chain2,mm2,icode2,seq2,nseq2)
      implicit none
      include 'pars.h'


      integer prealn(10,maxr),nprealn,prealn_len(10)   ! pre-defined alignment as initial seeds
      integer mm1(maxr),mm2(maxr),nseq1,nseq2
      character seq1(0:maxr),seq2(0:maxr)
      character chain1(maxr),chain2(maxr)
      character icode1(maxr),icode2(maxr)

      character*500 matchfile
      character*100 line, str
      character     AA1,AA2,tmpch
      character*5   matresid1(maxr),matresid2(maxr),nres1(maxr),nres2(maxr)
      character     matresch1(maxr),matresch2(maxr)

      integer i,ires,nmatch,n1,n2,match1,match2


      do i=1,nseq1
         if(icode1(i).eq." ")then
            write(nres1(i),'(I5)')mm1(i)
         else
            write(nres1(i),'(I4,A1)')mm1(i),icode1(i)
         endif
      end do
      do i=1,nseq2
         if(icode2(i).eq." ")then
            write(nres2(i),'(I5)')mm2(i)
         else
            write(nres2(i),'(I4,A1)')mm2(i),icode2(i)
         endif
      end do


      open(unit=16,file=matchfile,status='old')
      i=0
      nprealn=0
 5    read(16,1001,end=10)line

      if(line(1:6).eq.'PREALN')then
         nprealn=nprealn+1
         i=0
         nmatch=0
         do n1=1,nseq2
            prealn(nprealn,n1)=-1
         enddo
         goto 5
      endif
      if(line(1:3).eq.'END')then
         prealn_len(nprealn)=nmatch
         goto 5
      endif

      i=i+1
      read(line,'(I5)')ires

      if(ires.eq.i)then
         read(line,205)ires,matresch1(i),matresid1(i),AA1,tmpch,
     &        matresch2(i),matresid2(i),AA2
c         write(*,*)ires, matresch1(i), " ",matresid1(i)," ",AA1,
c     &        " - ", matresch2(i)," ",matresid2(i)," ",AA2

         match1=-1
         do n1=1,nseq1
            if(nres1(n1).eq.matresid1(i).and.chain1(n1).eq.matresch1(i)
     &           .and.seq1(n1).eq.AA1)then
               match1=n1
               goto 6
            end if
         end do
 6       continue

         match2=-1
         do n2=1,nseq2
            if(nres2(n2).eq.matresid2(i).and.chain2(n2).eq.matresch2(i)
     &           .and.seq2(n2).eq.AA2)then
               match2=n2
               goto 7
            end if
         end do
 7       continue

         if(match1.gt.0.and.match2.gt.0)then
            prealn(nprealn,match2)=match1
            nmatch=nmatch+1
            !print *,nprealn,nmatch,match1,match2
         endif

      endif
      goto 5
 10   continue

 205  format(i5,1x,a1,1x,a5,1x,a1,1x,a1,1x,a1,1x,a5,1x,a1)
 1001 format(A100)

      close(16)
      return
      end
ccccc==========================================================cccc

      subroutine readListFile(file,lst,num,npksel,pksel)
      use Constants
      implicit none

      character*500 file,lst(maxt),str1,str2
      character*2000 line
      character*30 pksel(maxt,maxp),tmpsel(maxp)
      character tab,tmp
      integer num,len,npksel(maxt)
      integer i,j,k

      tab=ACHAR(9)   !tabular key

      open(unit=17,file=file,status='old')
      num=0
 5    read(17,'(A)',end=10)line
      j=0
      k=0
      str1=''
      str2=''
      do i=1,len(line)
         if(line(i:i).ne.' '.and.line(i:i).ne.tab.and.k.eq.0)then
            j=j+1
            str1(j:j)=line(i:i)
         else if((line(i:i).eq.' '.or.line(i:i).eq.tab).and.j.gt.0.and.k.eq.0)then
            k=k+1
         else if(line(i:i).ne.' '.and.line(i:i).ne.tab.and.k.ge.1)then
            k=k+1
            str2(k:k)=line(i:i)
         else if(j.gt.0.and.k.gt.1) then
            exit
         endif
      enddo
      num=num+1
      lst(num)=trim(adjustl(str1))
      !print '(4A)','reading pocket selection: ',trim(adjustl(str1)),' : ',trim(adjustl(str2))

      if(k.gt.0)then
         call getPKind(adjustl(str2),tmpsel,npksel(num))
         do i=1,npksel(num)
            pksel(num,i)=tmpsel(i)
            !print '(2I4,1X,A)',num,i,trim(pksel(num,i))
         enddo
      else
         npksel(num)=0
      endif

      goto 5
 10   continue

      close(17)

      return
      end


ccccc==========================================================cccc
      subroutine mkPreAln(struct1,pkind1,struct2,pkind2,fullout,prealn)
      use Constants
      use PDBType
      use AlnResult
      implicit none

      type(pdbstruct) :: struct1,struct2
      type(alnres) :: fullout
      type(alignment) :: fullaln,prealn

      integer pkind1,pkind2
      integer i,j,k,naln,m1,m2, pkm1,pkm2
      integer map1(maxr),map2(maxr)

      fullaln = fullout%aln
      naln=0
      do i=1,fullaln%length
         m1=fullaln%map1(i)
         m2=fullaln%map2(i)

         pkm1=-1
         pkm2=-1
         do j=1,struct1%pk_len(pkind1)
            k=struct1%pk_res(pkind1,j)
            if(k.eq.m1)then
               pkm1=j
               exit
            endif
         enddo
         do j=1,struct2%pk_len(pkind2)
            k=struct2%pk_res(pkind2,j)
            if(k.eq.m2)then
               pkm2=j
               exit
            endif
         enddo

         if(pkm1.gt.0.and.pkm2.gt.0)then
            naln=naln+1
            map1(naln)=pkm1
            map2(naln)=pkm2
            !print *,naln,pkm1,pkm2
         endif
      enddo

      prealn%length=naln
      prealn%map1 = map1
      prealn%map2 = map2

      return
      end

ccccc===================================================ccccc
ccccc  P-value calculation for IS/TM-score.             ccccc
ccccc  Parameters were obtained empirically from one    ccccc
ccccc  million random interface alignments              ccccc
ccccc===================================================ccccc
      subroutine calcPvalue(score,nseq1,nseq2,pvalue,z)
      implicit none

      real*8  loc, scale, logt, logq
      real*8  z, pvalue
      real    score
      integer tlen, qlen, nseq1, nseq2

      qlen = min(nseq1,nseq2)   ! query/target is shorter
      tlen = max(nseq1,nseq2)   ! template is longer

      if(tlen.lt.5)tlen=5
      if(qlen.lt.5)qlen=5

      logt = log(dble(tlen))
      logq = log(dble(qlen))

      pvalue=-1
      z=-1


      !!! parameters based on Jeff's pockets
      loc   = 0.3117 + 0.0277*logt - 0.0290*logq
      scale = 0.0366 + 0.0025*logt - 0.0084*logq


      if(scale.lt.0.005) scale=0.005 !prevent insane values

      call calcEVDPV(score,loc,scale,pvalue,z)

      !write(*,'(F8.5,2I6,1x,F6.2,g12.4E3)')score,qlen,tlen,z,pvalue
      !write(*,'(2(A,F8.5))')'loc = ',loc,' scale = ',scale
      return
      end



ccccc===================================================ccccc
ccccc  P-value calculation for Gumbel Distribution.     ccccc
ccccc===================================================ccccc
      subroutine calcEVDPV(score,loc,scale,pvalue,z)

      real*8  loc, scale
      real*8  z, pvalue
      real    score

      z = (score - loc) / scale
      if( z .lt. 35 ) then
         pvalue = 1 -  dexp( - dexp( -z ) ) !Gumbel distribution
      else
c        below the double precision limit, use taylor expansion approaximation
         pvalue = dexp( -z )
      endif

      return
      end



ccccc====================================================cccccc
ccccc  read BLOSUM62 matrix
ccccc====================================================cccccc
      subroutine readBlosum62( seqmat )
      implicit none
      include 'pars.h'

      integer seqmat(0:19,0:19)

      seqmat(0,0:19)=(/6,0,0,-3,-3,-2,-4,-2,-3,-1,0,-4,-2,-2,-2,-2,-2,-3,-3,-2/) !G
      seqmat(1,0:19)=(/0,4,1,0,0,0,-1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-2,-2,-2,-3/) !A
      seqmat(2,0:19)=(/0,1,4,-1,-2,1,-2,-1,-1,0,1,-2,0,0,0,-1,-1,-2,-2,-3/) !S
      seqmat(3,0:19)=(/-3,0,-1,9,-1,-1,-1,-3,-1,-3,-3,-1,-3,-4,-3,-3,-3,-2,-2,-2/) !C
      seqmat(4,0:19)=(/-3,0,-2,-1,4,0,3,-2,1,-3,-3,1,-2,-2,-2,-3,-3,-1,-1,-3/) !V
      seqmat(5,0:19)=(/-2,0,1,-1,0,5,-1,-1,-1,-1,0,-1,-1,-1,-1,-1,-2,-2,-2,-2/) !T
      seqmat(6,0:19)=(/-4,-1,-2,-1,3,-1,4,-3,1,-3,-3,2,-3,-3,-3,-3,-3,0,-1,-3/) !I
      seqmat(7,0:19)=(/-2,-1,-1,-3,-2,-1,-3,7,-2,-1,-2,-3,-1,-1,-1,-2,-2,-4,-3,-4/) !P
      seqmat(8,0:19)=(/-3,-1,-1,-1,1,-1,1,-2,5,-3,-2,2,-1,-2,0,-1,-2,0,-1,-1/) !M
      seqmat(9,0:19)=(/-1,-2,0,-3,-3,-1,-3,-1,-3,6,1,-4,-1,2,0,-2,-1,-3,-3,-4/) !D
      seqmat(10,0:19)=(/0,-2,1,-3,-3,0,-3,-2,-2,1,6,-3,0,0,0,0,1,-3,-2,-4/) !N
      seqmat(11,0:19)=(/-4,-1,-2,-1,1,-1,2,-3,2,-4,-3,4,-2,-3,-2,-2,-3,0,-1,-2/) !L
      seqmat(12,0:19)=(/-2,-1,0,-3,-2,-1,-3,-1,-1,-1,0,-2,5,1,1,2,-1,-3,-2,-3/) !K
      seqmat(13,0:19)=(/-2,-1,0,-4,-2,-1,-3,-1,-2,2,0,-3,1,5,2,0,0,-3,-2,-3/) !E
      seqmat(14,0:19)=(/-2,-1,0,-3,-2,-1,-3,-1,0,0,0,-2,1,2,5,1,0,-3,-1,-2/) !Q
      seqmat(15,0:19)=(/-2,-1,-1,-3,-3,-1,-3,-2,-1,-2,0,-2,2,0,1,5,0,-3,-2,-3/) !R
      seqmat(16,0:19)=(/-2,-2,-1,-3,-3,-2,-3,-2,-2,-1,1,-3,-1,0,0,0,8,-1,2,-2/) !H
      seqmat(17,0:19)=(/-3,-2,-2,-2,-1,-2,0,-4,0,-3,-3,0,-3,-3,-3,-3,-1,6,3,1/) !F
      seqmat(18,0:19)=(/-3,-2,-2,-2,-1,-2,-1,-3,-1,-3,-2,-1,-2,-2,-1,-2,2,3,7,2/) !Y
      seqmat(19,0:19)=(/-2,-3,-3,-2,-3,-2,-3,-4,-1,-4,-4,-2,-3,-3,-2,-3,-2,1,2,11/) !W

      return
      end

ccccc====================================================cccccc
ccccc  uniform matrix for testing purpose
ccccc====================================================cccccc
      subroutine readUniMat( seqmat )
      implicit none
      include 'pars.h'

      integer seqmat(0:19,0:19),i,j

      do i=0,19
         do j=0,19
            seqmat(i,j)=1
         enddo
      enddo

      return
      end


ccccc====================================================cccccc
ccccc  reduced amino acid alphabet
ccccc  $str[0] = 'LVIMC';
cccc   $str[1] = 'AG';
cccc   $str[2] = 'ST';
cccc   $str[3] = 'P';
cccc   $str[4] = 'FYW';
cccc   $str[5] = 'EDNQ';
cccc   $str[6] = 'KR';
cccc   $str[7] = 'H';
ccccc====================================================cccccc
      subroutine readAA8( seqmat )
      implicit none
      include 'pars.h'

      integer seqmat(0:19,0:19)

      seqmat(0,0:19)=(/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/) !G
      seqmat(1,0:19)=(/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/) !A
      seqmat(2,0:19)=(/0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/) !S
      seqmat(3,0:19)=(/0,0,0,1,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0/) !C
      seqmat(4,0:19)=(/0,0,0,1,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0/) !V
      seqmat(5,0:19)=(/0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/) !T
      seqmat(6,0:19)=(/0,0,0,1,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0/) !I
      seqmat(7,0:19)=(/0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0/) !P
      seqmat(8,0:19)=(/0,0,0,1,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0/) !M
      seqmat(9,0:19)=(/0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0/) !D
      seqmat(10,0:19)=(/0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0/) !N
      seqmat(11,0:19)=(/0,0,0,1,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0/) !L
      seqmat(12,0:19)=(/0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0/) !K
      seqmat(13,0:19)=(/0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0/) !E
      seqmat(14,0:19)=(/0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0/) !Q
      seqmat(15,0:19)=(/0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0/) !R
      seqmat(16,0:19)=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0/) !H
      seqmat(17,0:19)=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1/) !F
      seqmat(18,0:19)=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1/) !Y
      seqmat(19,0:19)=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1/) !W


      return
      end

ccccc==========================================================cccc
ccccc==========================================================cccc
ccccc==========================================================cccc




ccccc==========================================================cccc
