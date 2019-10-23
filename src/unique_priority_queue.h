#include <set>
#include <queue>
#include <assert.h>

#ifndef UNIQUE_PRIORITY_QUEUE_H
#define UNIQUE_PRIORITY_QUEUE_H

using namespace std;

#define UNUSED(x) ((void)x)

const int    MAX_QUEUE_SIZE               = 6000000;
const int    MAX_NODES_MARGIN             = 500000;
const int    MAX_QUEUE_COPY_LENGTH        = 1000000;
const double QUEUE_COPY_LENGTH_PERCENTAGE = 1/((double)6);



template<class T>
struct do_nothing {
    void operator()(const T&){
        // intentionally left blank
    }
};

template<class T, class Container = vector<T>, class Compare = less<typename Container::value_type>>
class own_priority_queue : public priority_queue<T, Container, Compare> {
	public:
		Container& get_container() {
            return this->c;
        }
};

/**
 * Priority queue with unique (according to FuncCompare) elements of type T where the sorting is based on CostCompare.
 * If NDEBUG is *not* defined, there are some assertions that help catching errors in the provided comparision functions.
 */
template<class T, class CleanObsoleteElement = do_nothing<T>, class CostCompare = less<T>, class FuncCompare = CostCompare>
class unique_priority_queue
{
public:
    typedef typename own_priority_queue<T, vector<T>, CostCompare>::size_type size_type;

    /**
     * Return true if the element was inserted into the queue.
     * This happens if equivalent element is present or if the new element has a lower cost associated to it.
     * False is returned if no insertion into the queue took place.
     */
    bool push(const T& v)
    {
        const auto& insertion_pair = membership_.insert(v);
        if(insertion_pair.second)
        {
            queue_.push(v);
        }
        else if (CostCompare()(*(insertion_pair.first), v))
        {

            const auto number_erased = membership_.erase(*(insertion_pair.first));
            assert(number_erased == 1); UNUSED(number_erased);
            
            CleanObsoleteElement()(*(insertion_pair.first));
            const auto inserted = membership_.insert(v);
            assert(inserted.second); UNUSED(inserted);
            
            queue_ = own_priority_queue<T, vector<T>, CostCompare>();
            for(const auto& element : membership_) {
                queue_.push(element);
            }
            assert(queue_.size() == membership_.size());
            
            return true;
        }
        assert(queue_.size() == membership_.size());
        return insertion_pair.second;
    }

    bool push_delete_if_existing(const T& v) {
        bool success = push(v);
        if(!success) {
            CleanObsoleteElement()(v);
        }
        return success;
    }

    void pop()
    {
        assert(!queue_.empty() && queue_.size() == membership_.size());

        const auto& top_element   = queue_.top();
        const auto  number_erased = membership_.erase(top_element);

        assert(number_erased == 1); UNUSED(number_erased);

        queue_.pop();
        assert(queue_.size() == membership_.size());
    }

    const T& top() const
    {
        assert(!queue_.empty());
        return queue_.top();
    }

    bool empty() const
    {
        assert(queue_.size() == membership_.size());
        return queue_.empty();
    }

    size_type size() const
    {
        return queue_.size();
    }

    vector<T>& get_container() 
    {
        return queue_.get_container();
    }

    void delete_queue() 
    {
        vector<T> &v = get_container();
        for(typename vector<T>::const_iterator it = v.begin(); it != v.end();
            it++) {
            CleanObsoleteElement()(*it);
        }
        queue_      = own_priority_queue<T, vector<T>, CostCompare>();
        membership_ = set<T, FuncCompare>();
    }

    // clears the queue until a certain length is reached
    void update() 
    {	
        T temp_queue[MAX_QUEUE_SIZE];
        unsigned int length = min((int)(queue_.size() * QUEUE_COPY_LENGTH_PERCENTAGE), MAX_QUEUE_COPY_LENGTH);
        
        for(unsigned int i = 0; i < length; i++) {
            temp_queue[i] = queue_.top();
            queue_.pop();	

            last_node_copied = i;	
        }
        delete_queue();
        for(unsigned int i = 0; i < length; i++) {
            queue_.push(temp_queue[i]);		
        }
        
        cout << "RESULTING SIZE: " << queue_.size() << endl;
    }

    void restart(T& n) 
    {
        delete_queue();
        push(n);
    }
private:
    own_priority_queue<T, vector<T>, CostCompare> queue_;
    set<T, FuncCompare>                           membership_;
    unsigned int                                  last_node_copied = 0;
};
#endif
