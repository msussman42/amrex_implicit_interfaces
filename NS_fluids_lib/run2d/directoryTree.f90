! (c) Copyright Michael Metcalf and John Reid, 1999. This file may be
! freely used and copied for educational purposes provided this notice
! remains attached. Extracted from "Fortran 90/95 Explained" Oxford
! University Press (Oxford and New York), ISBN 0-19-851888-9.
!
!A recurring problem in computing is the need to manipulate a linked
!data structure.
!This might be a simple linked list like the one encountered in Section
!2.13, but often a more general tree structure is required.
!
!The example in this Appendix consists of a module that establishes and
!navigates one or more such trees, organized as a 'forest', and a short
!test program for it.  Here, each node is identified by a name and has
!any number of children, any number of siblings, and (optionally) some
!associated real data.  Each root node is regarded as having a common
!parent, the 'forest root' node, whose name is 'forest root'.  Thus,
!every node has a parent.  The module provides facilities for adding a
!named node to a specified parent, for enquiring about all the nodes
!that are offspring of a specified node, for removing a tree or subtree,
!and for performing I/O operations on a tree or subtree.
!
!The user-callable interfaces are:
!
!start:
!   must be called to initialize a forest.
!add_node:
!   stores the data provided at the node whose parent is specified
!   and sets up pointers to the parent and siblings (if any).
!remove_node:
!   deallocate all the storage occupied by a complete tree or
!   subtree.
!retrieve:
!   retrieves the data stored at a specified node and the names of
!   the parent and children.
!dump_tree:
!   write a complete tree or subtree.
!restore_tree:
!   read a complete tree or subtree.
!finish:
!   deallocate all the storage occupied by all the trees of the forest.
!
! changes by Anton VanderWyst
! University of Michigan
! Ann Arbor, MI 48015
!
! 07 Feb 2005 v2. changing 'char' node names to integer treenumbers.
!					Adding new structure fields
! 15 Feb 2005 v3. (v2 works). Changing to WhichPtInclude of position begin:end
!					instead of listing every point
! 01 Mar 2005 v4. (v3 works). Removing 'Charge' as a node:field
! 26 Apr 2005 v5. (v4 works). Removing unused variables, merging types

module directoryTree
!
! Strong typing imposed
	implicit none
!
! Only subroutine interfaces, the length of the character
! component, and the I/O unit number are public

	public  :: start, add_node, remove_node, retrieve,               &
		dump_tree, finish, node, find2, current, forest_root

! Module constants
	character(len=*), private, parameter :: eot = "End-of-Tree....."
	integer, parameter, public :: unit = 10,  & ! I/O unit number
	max_char = 16 ! length of character 
                                               ! component
! Define the basic tree type
	type node
		integer				:: name	! tree number
	! data for record
		real*8, pointer		:: ymin, ymax, xmin, xmax
		real*8, pointer		:: xcen, ycen
		real*8, pointer		:: aspectRatio
		integer, pointer	:: NumPtInclude, StartPosition, NumBoxes
		integer, dimension(:), pointer :: WhichChildren, WhichPtInclude

	! pointers to related records
		type(node), pointer			:: parent           ! parent node
		type(node), pointer			:: sibling          ! next sibling node
		type(node), pointer			:: child            ! first child node
	end type node

! Module variables
	type(node), pointer :: current     ! current node
	type(node), pointer :: forest_root ! the root of the forest
	integer, private    :: max_data    ! max size of data array
	integer, private, allocatable, target, dimension(:) :: names
                                           ! for returning list of names
! The module procedures

contains
	subroutine start ()
      ! Initialize the tree.
		allocate (forest_root)
		current => forest_root
		forest_root%name = 0
		nullify(forest_root%parent, forest_root%sibling, forest_root%child)
		allocate(forest_root%ymin); allocate(forest_root%ymax)
		allocate(forest_root%xmin); allocate(forest_root%xmax)
		allocate(forest_root%ycen); allocate(forest_root%xcen)
		allocate(forest_root%aspectRatio); allocate(forest_root%NumPtInclude)
		allocate(forest_root%StartPosition); allocate(forest_root%NumBoxes)
		allocate(forest_root%WhichPtInclude(0))
		allocate(forest_root%WhichChildren(0))

		max_data = 0
		allocate (names(0))
	end subroutine start

	subroutine find2(name)
		integer, intent(in) :: name
	  ! Make the module variable current point to the node with given name,
      ! or be null if the name is not there.
		type(node), pointer    :: root
      ! For efficiency, we search the tree rooted at current, and if this
      ! fails try its parent and so on until the forest root is reached.
		if (associated(current)) then
			root => current
			nullify(current)
		else
			root => forest_root
		endif

		do
			call look(root)
			if (associated(current)) then
				return
			end if
			root => root%parent
			if (.not.associated(root)) then
				exit
			end if
		end do
	contains 
		recursive subroutine look(root)
			type(node), pointer    :: root
		! (type(node), intent(in), target :: root   is standard conforming too)
		! Look for name in the tree rooted at root. If found, make the
		! module variable current point to the node
			type(node), pointer    :: child

			if (root%name == name) then
				current => root
			else
				child => root%child
				do
					if (.not.associated(child)) then
						exit
					end if
					call look(child)
					if (associated(current)) then
						return
					end if
					child => child%sibling
				end do
			end if
		end subroutine look
	end subroutine find2

	subroutine add_node(name, name_of_parent)
		integer, intent(in)             :: name, name_of_parent
	! For a root, name_of_parent = 0

	! Allocate a new tree node of type node, store the given name and
	! data there, set pointers to the parent and to its next sibling
	! (if any). If the parent is not found, the new node is treated as
	! a root. It is assumed that the node is not already present in the
	! forest.
		type(node), pointer :: new_node

		allocate (new_node)
		new_node%name = name
		allocate(new_node%ymin)
		allocate(new_node%ymax)
		allocate(new_node%xmin)
		allocate(new_node%xmax)
		allocate(new_node%xcen)
		allocate(new_node%ycen)
		allocate(new_node%aspectRatio)
		allocate(new_node%NumPtInclude)
		allocate(new_node%NumBoxes)
		allocate(new_node%StartPosition)
		allocate(new_node%WhichPtInclude(2))
		nullify(new_node%WhichChildren)

	! If name of parent is not null, search for it. 
	! If not found, print message.
		if (name_of_parent == 0) then
			current => forest_root
		else
			call find2 (name_of_parent)
			if (.not.associated(current)) then
				print *, "no parent ", name_of_parent, " found for ", name
				current => forest_root
			end if
		end if
		new_node%parent => current
		new_node%sibling => current%child
		current%child => new_node
		nullify(new_node%child)
	end subroutine add_node

	subroutine remove_node(name)
		integer, intent(in) :: name
	! Remove node and the subtree rooted on it (if any),
	! deallocating associated pointer targets.
		type(node), pointer :: parent, child, sibling

		call find2 (name)
		if (associated(current)) then
			parent =>  current%parent
			child => parent%child
			if (.not.associated(child, current)) then
		! Make it the first child, looping through the siblings to find it
		! and resetting the links
				parent%child => current
				sibling => child
				do
					if (associated (sibling%sibling, current)) then
						exit
					end if
					sibling => sibling%sibling
				end do
				sibling%sibling => current%sibling
				current%sibling => child
			end if
			call remove(current)
		end if
	end subroutine remove_node

	recursive subroutine remove (old_node)
	! Remove a first child node and the subtree rooted on it (if any),
	! deallocating associated pointer targets.
		type(node), pointer :: old_node
		type(node), pointer :: child, sibling

		child => old_node%child
		do
			if (.not.associated(child)) then
				exit
			end if
			sibling => child%sibling
			call remove(child)
			child => sibling
		end do
	! remove leaf node
		if (associated(old_node%parent)) then
			old_node%parent%child => old_node%sibling
		end if

		deallocate (old_node)
	end subroutine remove

	subroutine retrieve(name, data, parent, children)
		integer, intent(in)           :: name
		real, pointer, dimension(:)   :: data
		integer, intent(out)			:: parent
		integer, pointer, dimension(:) :: children
		! Returns a pointer to the data at the node, the name of the
		! parent, and a pointer to the names of the children.
		integer :: counter, i
		type(node), pointer :: child

		call find2 (name)
		if (associated(current)) then
         !data => current%y
			parent = current%parent%name
		! Count the number of children
			counter = 0
			child => current%child
			do
				if (.not.associated(child)) then
					exit
				end if
				counter = counter + 1
				child => child%sibling
			end do
			deallocate (names)
			allocate (names(counter))
		! and store their names
			children => names
			child => current%child
			do i = 1, counter
				children(i) = child%name
				child => child%sibling
			end do
		else
			nullify(data)
			parent = 0
			nullify(children)
		end if
	end subroutine retrieve

	subroutine dump_tree(root)
		integer, intent(in) :: root
	  ! Write out a complete tree followed by an end-of-tree record
	  ! unformatted on the file unit.
		call find2 (root)
		if (associated(current)) then
			call tree_out(current)
		end if
		write(unit = unit) eot, 0, eot
	contains 
		recursive subroutine tree_out (root)
		! Traverse a complete tree or subtree, writing out its contents
			type(node), intent(in) :: root     ! root node of tree
		! Local variable
			type(node), pointer    :: child

			write(unit = unit) root%name, root%parent%name
			child => root%child
			do
				if (.not.associated(child)) then
					exit
				end if
				call tree_out (child)
				child => child%sibling
			end do
		end subroutine tree_out
	end subroutine dump_tree

	subroutine finish ()
	! Deallocate all allocated targets.
		call remove (forest_root)
		deallocate(names)
	end subroutine finish

	SUBROUTINE change_node(ChangeNode, dataIn, dataChange, dataSlot)
		integer, intent(in)	:: dataChange
		real*8, intent(in)		:: dataIn
		integer, intent(in), optional :: dataSlot

	! Allocate a new tree node of type node, store the given name and
	! data there, set pointers to the parent and to its next sibling
	! (if any). CHANGE a particular data value in 'dataChange'.
	! If the parent is not found, the new node is treated as a root. 
		type(node), pointer		:: ChangeNode

		select case (dataChange)
    ! ymin=1, ymax=2, xmin=3, xmax=4, xcen=5, ycen=6, aspectRatio=7
	! NumPtInclude=8, WhichPtInclude=9, WhichChildren=10, 
	! StartPosition=11, NumBoxes=12
		case (1)
			ChangeNode%ymin=dataIn
		case (2)
			ChangeNode%ymax=dataIn
		case (3)
			ChangeNode%xmin=dataIn
		case (4)
			ChangeNode%xmax=dataIn
		case (5)
			ChangeNode%xcen=dataIn
		case (6)
			ChangeNode%ycen=dataIn
		case (7)
			ChangeNode%aspectRatio=dataIn
		case (8)
			ChangeNode%NumPtInclude=int(dataIn)
		case (9)
			ChangeNode%WhichPtInclude(dataSlot)=int(dataIn)
		case (10)
			ChangeNode%WhichChildren(dataSlot)=int(dataIn)
		case (11)
			ChangeNode%StartPosition=int(dataIn)
		case (12)
			ChangeNode%NumBoxes=int(dataIn)
		case default
			stop 'wrong "dataChange"'
		end select ! dataChange

	END SUBROUTINE change_node

end module directoryTree
